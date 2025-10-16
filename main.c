/*
 * Prime Finder via Fibonacci-Like Operator (C + GMP)
 *
 * Build: see Makefile
 *  - Enable OpenMP to parallelize heatmap/diagnostics and batched prime checks:
 *      CFLAGS += -fopenmp -flto
 *      LDFLAGS += -fopenmp -flto -lm -lgmp
 *  - Run with: OMP_NUM_THREADS=<n> ./prime_finder ...
 *
 * Notes:
 *  - Stage 1: same behavior, with early-abandon in window counting.
 *  - Heatmaps/diagnostics parallelized (OpenMP).
 *  - Stage 2: NEW batched + parallel primality testing. The FLO sequence
 *    is still generated sequentially, but candidate terms are collected
 *    into small chunks and tested concurrently across cores. Output
 *    (messages, stats) remains the same in spirit; progress still streams
 *    to stderr with rolling ETA based on observed per-check timing.
 */

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <stdint.h>
#include <inttypes.h>
#include <errno.h>
#ifdef _WIN32
#include <windows.h>
#include <direct.h>
#else
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef _WIN32
static long current_pid(void) { return (long)GetCurrentProcessId(); }
#else
static long current_pid(void) { return (long)getpid(); }
#endif

/* ---------- Expected-trials estimator (unchanged API) ---------- */

static inline double _prime_density_correction(double lnN)
{
    const double inv = 1.0 / lnN;
    const double inv2 = inv * inv;
    const double inv3 = inv2 * inv;
    return 1.0 + inv + 2.0 * inv2 + 6.0 * inv3;
}

/* Expected primality tests to get one prime with ~digits digits. */
double expected_prime_trials(int digits, int odd_only, int use_asymptotic_corrections)
{
    if (digits <= 0)
        return NAN;
    const double LN10 = 2.30258509299404568402;
    const double lnN = digits * LN10;

    double density = 1.0 / lnN;
    if (use_asymptotic_corrections)
        density *= _prime_density_correction(lnN);
    if (odd_only)
        density *= 2.0;

    const double expected = 1.0 / density;
    if (!isfinite(expected) || expected <= 0.0)
        return NAN;
    return expected;
}

/* ---------- Types & helpers ---------- */

static void configure_streams_for_pm2(int argc, char **argv)
{
    /* Ensure logs reach PM2 immediately even under non-tty buffering. */
    setvbuf(stdout, NULL, _IOLBF, 0);
    setvbuf(stderr, NULL, _IONBF, 0);

    fprintf(stderr, "[PM2] stdout=line-buffered stderr=unbuffered (pid=%ld)\n", current_pid());
    fprintf(stderr, "[PM2] argv:");
    for (int i = 0; i < argc; ++i)
        fprintf(stderr, " %s", argv[i]);
    fprintf(stderr, "\n");
}

static double wall_time_now(void)
{
#ifdef _WIN32
    static LARGE_INTEGER freq;
    static int freq_initialized = 0;
    LARGE_INTEGER counter;
    if (!freq_initialized)
    {
        QueryPerformanceFrequency(&freq);
        freq_initialized = 1;
    }
    QueryPerformanceCounter(&counter);
    return (double)counter.QuadPart / (double)freq.QuadPart;
#else
#ifdef CLOCK_MONOTONIC
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
#else
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec + (double)tv.tv_usec * 1e-6;
#endif
#endif
}

typedef struct
{
    int s1, s2;
    int hits;
    double density;
} SeedStat;

static const unsigned int SMALL_PRIMES[] = {3, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59};
static const size_t SMALL_PRIME_COUNT = sizeof(SMALL_PRIMES) / sizeof(SMALL_PRIMES[0]);

static int passes_small_prime_filters(const mpz_t n)
{
    if (mpz_sgn(n) <= 0)
        return 0;
    if (mpz_even_p(n))
        return 0;
    if (mpz_tdiv_ui(n, 5) == 0)
        return 0;
    for (size_t i = 0; i < SMALL_PRIME_COUNT; ++i)
    {
        if (mpz_tdiv_ui(n, SMALL_PRIMES[i]) == 0)
            return 0;
    }
    return 1;
}

static void flo_next(mpz_t a, mpz_t b)
{
    // (a,b) -> (b, a+b)
    mpz_t tmp;
    mpz_init(tmp);
    mpz_add(tmp, a, b);
    mpz_set(a, b);
    mpz_set(b, tmp);
    mpz_clear(tmp);
}

/* count hits with early-abandon vs current best (threshold in [0,1]) */
static int count_primes_in_window_abandon(int s1, int s2, int window, int best_hits_so_far, double threshold, int *abandoned_out)
{
    mpz_t a, b;
    mpz_init_set_ui(a, (unsigned long)s1);
    mpz_init_set_ui(b, (unsigned long)s2);
    int hits = 0;
    int abandoned = 0;

    /* term 0 */
    if (!(s1 > 10 && (s1 % 2 == 0 || s1 % 5 == 0)))
    {
        if (mpz_probab_prime_p(a, 10) > 0)
            hits++;
    }
    if (best_hits_so_far > 0)
    {
        int remaining = window - 1;
        if (hits + remaining < (int)ceil(threshold * (double)best_hits_so_far))
        {
            abandoned = 1;
            goto FINISH;
        }
    }

    /* term 1 */
    if (!(s2 > 10 && (s2 % 2 == 0 || s2 % 5 == 0)))
    {
        if (mpz_probab_prime_p(b, 10) > 0)
            hits++;
    }
    if (best_hits_so_far > 0)
    {
        int remaining = window - 2;
        if (hits + remaining < (int)ceil(threshold * (double)best_hits_so_far))
        {
            abandoned = 1;
            goto FINISH;
        }
    }

    for (int i = 2; i < window; i++)
    {
        flo_next(a, b);
        if (!(mpz_divisible_ui_p(b, 2) || mpz_divisible_ui_p(b, 5)))
        {
            if (mpz_probab_prime_p(b, 10) > 0)
                hits++;
        }
        if (best_hits_so_far > 0)
        {
            int remaining = window - (i + 1);
            if (hits + remaining < (int)ceil(threshold * (double)best_hits_so_far))
            {
                abandoned = 1;
                break;
            }
        }
    }

FINISH:
    mpz_clear(a);
    mpz_clear(b);
    if (abandoned_out)
        *abandoned_out = abandoned;
    return hits;
}

static int cmp_seedstat(const void *A, const void *B)
{
    const SeedStat *a = (const SeedStat *)A;
    const SeedStat *b = (const SeedStat *)B;
    if (a->hits != b->hits)
        return (b->hits - a->hits);
    if (a->density < b->density)
        return 1;
    if (a->density > b->density)
        return -1;
    return 0;
}

/* ---------- Minimal PPM heatmap writer (no extra deps) ---------- */
/* We generate PPM (P6) images representing prime-hit diagnostics:
   - heatmap.ppm             : standard prime hits
   - heatmap_gcd.ppm         : gcd(s1,s2) intensity map (g>1 bands)
   - heatmap_exclusions.ppm  : fraction of quick-excluded terms (mod 2 or 5) in the window
   X-axis: seed2, Y-axis: seed1.                                      */

static void palette_fire(double t, unsigned char *r, unsigned char *g, unsigned char *b)
/* t in [0,1]; maps black→red→orange→yellow→white */
{
    if (t < 0)
        t = 0;
    if (t > 1)
        t = 1;
    double x = t;
    double R, G, B;
    if (x < 0.25)
    {
        R = 4.0 * x;
        G = 0.0;
        B = 0.0;
    }
    else if (x < 0.50)
    {
        R = 1.0;
        G = 4.0 * (x - 0.25);
        B = 0.0;
    }
    else if (x < 0.75)
    {
        R = 1.0;
        G = 1.0;
        B = 0.0;
    }
    else
    {
        R = 1.0;
        G = 1.0;
        B = 4.0 * (x - 0.75);
        if (B > 1.0)
            B = 1.0;
    }
    *r = (unsigned char)(255.0 * R);
    *g = (unsigned char)(255.0 * G);
    *b = (unsigned char)(255.0 * B);
}

static void write_heatmap_ppm(const SeedStat *grid_stats,
                              int seed_min, int seed_max,
                              const char *path)
{
    int W = (seed_max - seed_min + 1);
    int H = (seed_max - seed_min + 1);
    int N = W * H;

    /* parallel min/max reduction */
    int min_hits = grid_stats[0].hits, max_hits = grid_stats[0].hits;
#pragma omp parallel for reduction(min : min_hits) reduction(max : max_hits) schedule(static) if (N > 2048)
    for (int i = 0; i < N; ++i)
    {
        int h = grid_stats[i].hits;
        if (h < min_hits)
            min_hits = h;
        if (h > max_hits)
            max_hits = h;
    }
    double span = (double)(max_hits - min_hits);
    if (span <= 0)
        span = 1.0;

    FILE *fp = fopen(path, "wb");
    if (!fp)
    {
        fprintf(stderr, "[HEATMAP] Failed to open %s for writing\n", path);
        return;
    }

    fprintf(fp, "P6\n%d %d\n255\n", W, H);

    /* Write rows in order; parallelize within row to keep deterministic file layout */
    for (int row = 0; row < H; ++row)
    {
        unsigned char *rowbuf = (unsigned char *)malloc(3 * W);
#pragma omp parallel for schedule(static) if (W > 64)
        for (int col = 0; col < W; ++col)
        {
            int idx = row * W + col;
            double t = ((double)grid_stats[idx].hits - (double)min_hits) / span;
            unsigned char r, g, b;
            palette_fire(t, &r, &g, &b);
            rowbuf[3 * col + 0] = r;
            rowbuf[3 * col + 1] = g;
            rowbuf[3 * col + 2] = b;
        }
        fwrite(rowbuf, 1, 3 * W, fp);
        free(rowbuf);
    }
    fclose(fp);

    /* CSV */
    FILE *csv = fopen("heatmap.csv", "w");
    if (csv)
    {
        fprintf(csv, "# rows = seed1 [%d..%d], cols = seed2 [%d..%d], values = prime hits in first 'window' terms\n",
                seed_min, seed_max, seed_min, seed_max);
        for (int row = 0; row < H; ++row)
        {
            for (int col = 0; col < W; ++col)
            {
                int idx = row * W + col;
                fprintf(csv, "%d%s", grid_stats[idx].hits, (col == W - 1) ? "" : ",");
            }
            fprintf(csv, "\n");
        }
        fclose(csv);
    }

    printf("[HEATMAP] Wrote %s (P6) and heatmap.csv\n", path);
}

/* gcd heatmap: brightness increases with gcd(s1,s2); 0 if gcd==1 */
static void write_gcd_heatmap_ppm(int seed_min, int seed_max, const char *path)
{
    int W = (seed_max - seed_min + 1);
    int H = (seed_max - seed_min + 1);
    double max_val = log2((double)seed_max);

    FILE *fp = fopen(path, "wb");
    if (!fp)
    {
        fprintf(stderr, "[HEATMAP] Failed to open %s for writing\n", path);
        return;
    }
    fprintf(fp, "P6\n%d %d\n255\n", W, H);

    for (int row = 0; row < H; ++row)
    {
        unsigned char *rowbuf = (unsigned char *)malloc(3 * W);
        int s1 = seed_min + row;
#pragma omp parallel for schedule(static) if (W > 64)
        for (int col = 0; col < W; ++col)
        {
            int s2 = seed_min + col;
            int a = s1, b = s2;
            while (b)
            {
                int t = a % b;
                a = b;
                b = t;
            }
            int g = a;

            double t = 0.0;
            if (g > 1 && max_val > 0.0)
                t = log2((double)g) / max_val;
            unsigned char r, gc, bc;
            palette_fire(t, &r, &gc, &bc);

            rowbuf[3 * col + 0] = r;
            rowbuf[3 * col + 1] = gc;
            rowbuf[3 * col + 2] = bc;
        }
        fwrite(rowbuf, 1, 3 * W, fp);
        free(rowbuf);
    }
    fclose(fp);
    printf("[HEATMAP] Wrote %s (gcd intensity)\n", path);
}

/* exclusions heatmap: fraction of terms (in first 'window') eliminated by quick screens (mod 2 or 5) */
static void write_exclusions_heatmap_ppm(int seed_min, int seed_max, int window, const char *path)
{
    int W = (seed_max - seed_min + 1);
    int H = (seed_max - seed_min + 1);

    FILE *fp = fopen(path, "wb");
    if (!fp)
    {
        fprintf(stderr, "[HEATMAP] Failed to open %s for writing\n", path);
        return;
    }
    fprintf(fp, "P6\n%d %d\n255\n", W, H);

    for (int row = 0; row < H; ++row)
    {
        unsigned char *rowbuf = (unsigned char *)malloc(3 * W);
        int s1 = seed_min + row;

#pragma omp parallel for schedule(static) if (W > 64)
        for (int col = 0; col < W; ++col)
        {
            int s2 = seed_min + col;

            mpz_t a, b;
            mpz_init_set_ui(a, (unsigned long)s1);
            mpz_init_set_ui(b, (unsigned long)s2);
            int excl = 0, terms = 0;

            terms++;
            if ((s1 > 10) && (s1 % 2 == 0 || s1 % 5 == 0))
                excl++;
            terms++;
            if ((s2 > 10) && (s2 % 2 == 0 || s2 % 5 == 0))
                excl++;

            for (int i = 2; i < window; ++i)
            {
                flo_next(a, b);
                terms++;
                if (mpz_divisible_ui_p(b, 2) || mpz_divisible_ui_p(b, 5))
                    excl++;
            }

            mpz_clear(a);
            mpz_clear(b);

            double frac = (terms > 0) ? ((double)excl / (double)terms) : 0.0;
            unsigned char r, g, bc;
            palette_fire(frac, &r, &g, &bc);
            rowbuf[3 * col + 0] = r;
            rowbuf[3 * col + 1] = g;
            rowbuf[3 * col + 2] = bc;
        }
        fwrite(rowbuf, 1, 3 * W, fp);
        free(rowbuf);
    }
    fclose(fp);
    printf("[HEATMAP] Wrote %s (quick-exclusion fraction)\n", path);
}

/* ---------- Stage 2 parallel check support ---------- */

typedef struct
{
    mpz_t n;
    int term_index; /* term counter within the post-threshold loop */
} Candidate;

#ifndef PAR_CHUNK_SIZE
#define PAR_CHUNK_SIZE 1 /* number of candidate terms per parallel chunk */
#endif

#ifndef MR_REPS_CAND
#define MR_REPS_CAND 8
#endif

#ifndef MR_REPS_FINAL
#define MR_REPS_FINAL 25
#endif



/* Fill up to PAR_CHUNK_SIZE candidates by advancing (a,b) terms.
   Applies small-prime filters before enqueueing. Returns number of
   candidates filled, and advances t (term counter). */
static int build_candidate_chunk(Candidate *cand, mpz_t a, mpz_t b, int *t, int max_terms)
{
    int m = 0;
    while (m < PAR_CHUNK_SIZE && (max_terms < 0 || *t < max_terms))
    {
        /* current term is in b; quick exclusions */
        if (passes_small_prime_filters(b))
        {
            mpz_init_set(cand[m].n, b);
            cand[m].term_index = *t;
            m++;
        }
        /* advance to next term */
        flo_next(a, b);
        (*t)++;
    }
    return m;
}

/* Clear a candidate chunk */
static void clear_candidate_chunk(Candidate *cand, int m)
{
    for (int i = 0; i < m; ++i)
    {
        mpz_clear(cand[i].n);
    }
}

/* Parallel primality test over a chunk. Returns:
   - found_idx >=0 : index within chunk of first prime found (notionally earliest term)
   - else -1       : none prime in this chunk.
   We ensure "earliest term" semantics by resolving ties after parallel section. */
static int parallel_test_chunk(Candidate *cand, int m)
{
    if (m <= 0)
        return -1;

    int *is_prime = (int *)calloc(m, sizeof(int));
    if (!is_prime)
        return -1;

#pragma omp parallel for schedule(static) if (m > 1)
    for (int i = 0; i < m; ++i)
    {
        int pr = mpz_probab_prime_p(cand[i].n, MR_REPS_CAND);
        is_prime[i] = (pr > 0) ? 1 : 0;
    }

    /* Find the smallest term_index that is prime (earliest in sequence) */
    int found_idx = -1;
    int best_term = 0x7fffffff;
    for (int i = 0; i < m; ++i)
    {
        if (is_prime[i])
        {
            if (cand[i].term_index < best_term)
            {
                best_term = cand[i].term_index;
                found_idx = i;
            }
        }
    }

    free(is_prime);
    return found_idx;
}

/* ---------- Main ---------- */

int main(int argc, char **argv)
{
    configure_streams_for_pm2(argc, argv);

    double program_t0 = wall_time_now();

    /* CLI-compatible defaults (unchanged) */
    int seed_min = 0, seed_max = 0;
    const int window = 500;
    int top = 10, max_terms = -1, target_digits = 10000;
    const char *out_path = "results.txt";
    const char *batch_path = "batch-results.txt";

    for (int i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "--target-digits") && i + 1 < argc)
            target_digits = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--out") && i + 1 < argc)
            out_path = argv[++i];
        else if (!strcmp(argv[i], "--batch") && i + 1 < argc)
            batch_path = argv[++i];
    }

    unsigned int entropy = (unsigned int)time(NULL);
    entropy ^= (unsigned int)current_pid();
    entropy ^= (unsigned int)(uintptr_t)&entropy;
    srand(entropy);
    const int seed_span = 10000 - 1;
    seed_min = 1 + (rand() % seed_span);
    seed_max = seed_min + 200;

    int omp_threads = 1;
#ifdef _OPENMP
#ifdef DEFAULT_OMP_THREADS
    omp_set_num_threads(DEFAULT_OMP_THREADS);
#endif
#pragma omp parallel
    {
#pragma omp master
        omp_threads = omp_get_num_threads();
    }
    fprintf(stderr, "[OMP] Using up to %d threads for heatmap/diagnostics and parallel prime checks\n", omp_threads);
#endif

    printf("[RANDOM] seed_min=%d seed_max=%d (span=%d) window=%d\n",
           seed_min, seed_max, seed_max - seed_min + 1, window);

    char max_terms_label[32];
    if (max_terms < 0)
        strcpy(max_terms_label, "unbounded");
    else
        snprintf(max_terms_label, sizeof(max_terms_label), "%d", max_terms);

    printf("== FLO Prime Finder (C + GMP) ==\n");
    printf("seed_min=%d seed_max=%d window=%d top=%d target_digits=%d max_terms=%s\n",
           seed_min, seed_max, window, top, target_digits, max_terms_label);

    /* Stage 1: Seed scan with optional on-disk cache */
    int total = (seed_max - seed_min + 1) * (seed_max - seed_min + 1);
    SeedStat *stats = (SeedStat *)malloc(sizeof(SeedStat) * total);
    SeedStat *stats_grid = (SeedStat *)malloc(sizeof(SeedStat) * total); /* preserve pre-sort order for heatmap */

    double seed_scan_seconds = 0.0;
    const double keepup_threshold = 0.80; /* abandon if cannot reach 80% of current best */
    int best_hits_so_far = 0;
    int abandoned_count = 0;
    int idx = 0;
    double next_progress = 0.05;

    double scan_t0 = wall_time_now();
    printf("[STATS] target: %d | scanning seeds: 0%% complete\n", target_digits);
    for (int s1 = seed_min; s1 <= seed_max; ++s1)
    {
        for (int s2 = seed_min; s2 <= seed_max; ++s2)
        {
            int abandoned_flag = 0;
            int hits = count_primes_in_window_abandon(s1, s2, window, best_hits_so_far, keepup_threshold, &abandoned_flag);

            stats[idx].s1 = s1;
            stats[idx].s2 = s2;
            stats[idx].hits = hits;
            stats[idx].density = (double)hits / (double)(window > 0 ? window : 1);

            stats_grid[idx] = stats[idx];
            idx++;

            double progress = (double)idx / (double)total;
            while (progress >= next_progress && next_progress < 1.001)
            {
                int pct = (int)(next_progress * 100.0 + 0.5);
                printf("[STATS] target: %d | scanning seeds: %d%% complete\n", target_digits, pct);
                next_progress += 0.05;
            }

            if (abandoned_flag)
                abandoned_count++;
            if (hits > best_hits_so_far)
                best_hits_so_far = hits;
        }
    }
    if (next_progress <= 1.001)
        printf("[STATS] target: %d | scanning seeds: 100%% complete\n", target_digits);
    seed_scan_seconds = wall_time_now() - scan_t0;
    printf("[STATS] scanned %d seed pairs in %.2fs | early-abandoned=%d (threshold=%.2f of best)\n",
           total, seed_scan_seconds, abandoned_count, keepup_threshold);

    /* Heatmaps/diagnostics */
    write_heatmap_ppm(stats_grid, seed_min, seed_max, "heatmap.ppm");
    write_gcd_heatmap_ppm(seed_min, seed_max, "heatmap_gcd.ppm");
    write_exclusions_heatmap_ppm(seed_min, seed_max, window, "heatmap_exclusions.ppm");

    int gcd_gt1 = 0, diagonal = 0;
    for (int i = 0; i < total; ++i)
    {
        int a = stats_grid[i].s1, b = stats_grid[i].s2;
        int x = a, y = b;
        while (y)
        {
            int t = x % y;
            x = y;
            y = t;
        }
        int g = x;
        if (g > 1)
            gcd_gt1++;
        if (stats_grid[i].s1 == stats_grid[i].s2)
            diagonal++;
    }
    printf("[DIAG] gcd(s1,s2)>1 seeds: %d/%d | diagonal seeds: %d\n", gcd_gt1, total, diagonal);

    /* Leaderboard */
    qsort(stats, total, sizeof(SeedStat), cmp_seedstat);

    FILE *batch = fopen(batch_path, "w");
    fprintf(batch, "# Seed Scan Leaderboard\n# window=%d\n", window);
    fprintf(batch, "rank,s1,s2,hits,density\n");
    for (int i = 0; i < (top * 2 < total ? top * 2 : total); ++i)
    {
        fprintf(batch, "%d,%d,%d,%d,%.3f\n", i + 1, stats[i].s1, stats[i].s2, stats[i].hits, stats[i].density);
        printf("%d,%d,%d,%d,%.3f\n", i + 1, stats[i].s1, stats[i].s2, stats[i].hits, stats[i].density);
    }
    fclose(batch);

    /* choose seed pair with MAX hits; randomize among ties */
    srand((unsigned)time(NULL));
    int max_hits = stats[0].hits; /* after sort, stats[0] has the maximum hits */
    int tie_count = 0;
    while (tie_count < total && stats[tie_count].hits == max_hits)
    {
        tie_count++;
    }
    int chosen_idx = (tie_count > 0) ? (rand() % tie_count) : 0;

    if (chosen_idx != 0)
    {
        SeedStat tmp = stats[0];
        stats[0] = stats[chosen_idx];
        stats[chosen_idx] = tmp;
    }
    int old_top = top;
    top = 1; /* Force single-seed focused search */
    printf("[SELECT] Using seed with max prime hits (ties randomized): (%d,%d) hits=%d (tie-group=%d)\n",
           stats[0].s1, stats[0].s2, stats[0].hits, tie_count);

    /* Stage 2: Focused search with PARALLEL primality testing */
    FILE *out = fopen(out_path, "w");
    fprintf(out, "# Focused FLO prime search\n# target_digits=%d max_terms=%s\n", target_digits, max_terms_label);

    printf("[TIMING] seed_scan_elapsed=%.3fs\n", seed_scan_seconds);
    fprintf(out, "# seed_scan_elapsed=%.3fs\n", seed_scan_seconds);

    double expected_checks = expected_prime_trials(target_digits, /*odd_only=*/1, /*corrections=*/0);
    printf("Estimated primality tests per prime @ %d digits (odd-only): %.1f\n", target_digits, expected_checks);
    fprintf(out, "Estimated primality tests per prime @ %d digits (odd-only): %.1f\n", target_digits, expected_checks);

    mpz_t a, b;
   mpz_init(a);
   mpz_init(b);

   int total_found = 0;
    double total_search_seconds = 0.0;

    for (int i = 0; i < (top < total ? top : total); ++i)
    {
        int s1 = stats[i].s1, s2 = stats[i].s2;
        fprintf(out, "\n== Seed (%d,%d) rank=%d window_hits=%d density=%.3f\n",
                s1, s2, i + 1, stats[i].hits, stats[i].density);
        printf("[SEARCH] Seed (%d,%d)\n", s1, s2);

        double seed_search_start = wall_time_now();
        double seed_search_elapsed = 0.0;

        /* Initialize sequence */
        mpz_set_ui(a, (unsigned long)s1);
        mpz_set_ui(b, (unsigned long)s2);

        /* warm-up to target digits (no prime checks before threshold) */
        int terms_advanced = 0;
        while (mpz_sizeinbase(b, 10) < (unsigned)target_digits)
        {
            flo_next(a, b);
            terms_advanced++;
        }
        
        printf("Advanced %d terms to reach %d-digit threshold (current digits=%lu)\n",
               terms_advanced, target_digits, mpz_sizeinbase(b, 10));

        double search_phase_start = wall_time_now();
        long long checks = 0;
        int found_for_seed = 0;
/* NEW: batched + parallel prime checking */
        int t = 0; /* term counter post-threshold */
        while (1)
        {
            Candidate cand[PAR_CHUNK_SIZE];
            int m = build_candidate_chunk(cand, a, b, &t, max_terms);
            if (m <= 0)
                continue;

            double chunk_t0 = wall_time_now();
            int found_idx = parallel_test_chunk(cand, m);
            double chunk_t1 = wall_time_now();

            /* Update checks (exact): if found, only count up to the winning candidate */
            if (found_idx >= 0)
                checks += (found_idx + 1);
            else
                checks += m;

            /* Derive per-check timing + ETA from chunk stats */
            double chunk_secs = chunk_t1 - chunk_t0;
            double elapsed_secs = chunk_t1 - search_phase_start;
            double avg_chk = (checks > 0) ? (elapsed_secs / (double)checks) : (m > 0 ? chunk_secs / (double)m : 0.0);
            double remain = (expected_checks > (double)checks) ? (expected_checks - (double)checks) : 0.0;
            double eta_sec = remain * avg_chk;
            double eta_parallel = eta_sec;
            if (omp_threads > 1)
                eta_parallel = eta_sec / (double)omp_threads;

            fprintf(stderr,
                    "\rcheck #%lld | last≈%.6fs avg=%.6fs | est. remain≈%.1f checks | ETA≈%.1fs (@%d thr)   ",
                    checks,
                    (m > 0 ? (chunk_secs / (double)m) : 0.0),
                    avg_chk,
                    remain,
                    eta_parallel,
                    omp_threads);

            if (found_idx >= 0)
            {
                char *prime_str = mpz_get_str(NULL, 10, cand[found_idx].n);
                if (!prime_str)
                {
                    clear_candidate_chunk(cand, m);
                    continue;
                }

                size_t str_len = strlen(prime_str);
                size_t clip = 12;
                char preview[64];
                if (str_len > clip * 2 + 4)
                {
                    char head[clip + 1];
                    char tail[clip + 1];
                    memcpy(head, prime_str, clip);
                    head[clip] = '\0';
                    memcpy(tail, prime_str + (str_len - clip), clip);
                    tail[clip] = '\0';
                    snprintf(preview, sizeof(preview), "%s...%s", head, tail);
                }
                else
                {
                    snprintf(preview, sizeof(preview), "%s", prime_str);
                }

                printf("[CANDIDATE] seed=(%d,%d) digits=%zu idx=%d checks=%lld preview=%s\n",
                       s1, s2, str_len, cand[found_idx].term_index, checks, preview);

                int final_pr = mpz_probab_prime_p(cand[found_idx].n, MR_REPS_FINAL);
                if (final_pr <= 0)
                {
                    printf("[CANDIDATE] Rejected after MR_REPS_FINAL=%d rounds (composite detected)\n", MR_REPS_FINAL);
                    free(prime_str);
                    clear_candidate_chunk(cand, m);
                    continue;
                }

                /* stop timer on confirmed prime */
                double secs = wall_time_now() - search_phase_start;
                int dcount = (int)str_len;
                double eff_ratio = (checks > 0) ? (expected_checks / (double)checks) : 0.0;
                double expected_checks_found = expected_prime_trials(dcount, /*odd_only=*/1, /*corrections=*/0);
                double eff_ratio_found = (checks > 0) ? (expected_checks_found / (double)checks) : 0.0;

                printf("FOUND probable prime | seed=(%d,%d) digits=%d\n", s1, s2, dcount);
                printf("  time_to_find: %.3f sec | checks: %lld | expected_checks: %.1f | efficiency (expected/actual): %.3f\n",
                       secs, checks, expected_checks, eff_ratio);
                printf("  expected_checks_actual_digits(%d): %.1f | efficiency_actual (expected/actual): %.3f\n",
                       dcount, expected_checks_found, eff_ratio_found);
                printf("  prime: %s\n", prime_str);

                fprintf(out, "FOUND probable prime | seed=(%d,%d) digits=%d\n", s1, s2, dcount);
                fprintf(out, "  time_to_find: %.3f sec | checks: %lld | expected_checks: %.1f | efficiency (expected/actual): %.3f\n",
                        secs, checks, expected_checks, eff_ratio);
                fprintf(out, "  expected_checks_actual_digits(%d): %.1f | efficiency_actual (expected/actual): %.3f\n",
                        dcount, expected_checks_found, eff_ratio_found);
                fprintf(out, "  prime: %s\n", prime_str);
                fflush(out);

                free(prime_str);
                clear_candidate_chunk(cand, m);
                total_found++;
                found_for_seed = 1;
                break;
            }

            /* no confirmed prime in this chunk; clear and continue */
            clear_candidate_chunk(cand, m);
        }

        if (!found_for_seed && max_terms >= 0)
        {
            printf("  No prime found within %d checks-window (post-threshold) for this seed.\n", max_terms);
            fprintf(out, "  No prime found within %d checks-window (post-threshold).\n", max_terms);
        }

        seed_search_elapsed = wall_time_now() - seed_search_start;
        total_search_seconds += seed_search_elapsed;
        printf("[TIMING] seed=(%d,%d) search_elapsed=%.3fs\n", s1, s2, seed_search_elapsed);
        fprintf(out, "# seed=(%d,%d) search_elapsed=%.3fs\n", s1, s2, seed_search_elapsed);
    }

    double program_elapsed = wall_time_now() - program_t0;
    fprintf(out, "# timing_totals seed_scan=%.3fs prime_search=%.3fs overall=%.3fs\n",
            seed_scan_seconds, total_search_seconds, program_elapsed);
    printf("[TIMING] totals: seed_scan=%.3fs prime_search=%.3fs overall=%.3fs\n",
           seed_scan_seconds, total_search_seconds, program_elapsed);

    fprintf(out, "\n# total_found=%d\n", total_found);
    fclose(out);

    mpz_clear(a);
    mpz_clear(b);
    free(stats_grid);
    free(stats);

    (void)old_top; /* preserve prior variable, no-op */

    return 0;
}
