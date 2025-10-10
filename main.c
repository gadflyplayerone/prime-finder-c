/*
 * Prime Finder via Fibonacci-Like Operator (C + GMP)
 *
 * Build: see Makefile
 */

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

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

typedef struct
{
    int s1, s2;
    int hits;
    double density;
} SeedStat;

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
    /* possible early check after term 0 */
    if (best_hits_so_far > 0)
    {
        int remaining = window - 1; /* terms left to observe */
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

static int count_primes_in_window(int s1, int s2, int window)
{
    return count_primes_in_window_abandon(s1, s2, window, /*best_hits_so_far=*/0, /*threshold=*/0.0, /*abandoned_out=*/NULL);
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
   - heatmap.ppm           : standard prime hits
   - heatmap_gcd.ppm      : gcd(s1,s2) intensity map (g>1 bands)
   - heatmap_exclusions.ppm: fraction of quick-excluded terms (mod 2 or 5) in the window
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

    int min_hits = grid_stats[0].hits, max_hits = grid_stats[0].hits;
    for (int i = 0; i < W * H; ++i)
    {
        if (grid_stats[i].hits < min_hits)
            min_hits = grid_stats[i].hits;
        if (grid_stats[i].hits > max_hits)
            max_hits = grid_stats[i].hits;
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
    for (int row = 0; row < H; ++row)
    {
        for (int col = 0; col < W; ++col)
        {
            int idx = row * W + col;
            double t = ((double)grid_stats[idx].hits - (double)min_hits) / span;
            unsigned char r, g, b;
            palette_fire(t, &r, &g, &b);
            unsigned char px[3] = {r, g, b};
            fwrite(px, 1, 3, fp);
        }
    }
    fclose(fp);

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

    /* find max gcd to normalize (upper bound is seed_max) */
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
        int s1 = seed_min + row;
        for (int col = 0; col < W; ++col)
        {
            int s2 = seed_min + col;
            int g = 1;
            /* fast gcd */
            int a = s1, b = s2;
            while (b)
            {
                int t = a % b;
                a = b;
                b = t;
            }
            g = a;

            double t = 0.0;
            if (g > 1 && max_val > 0.0)
                t = log2((double)g) / max_val; /* [0,1] */
            unsigned char r, gC, bC;
            palette_fire(t, &r, &gC, &bC);
            unsigned char px[3] = {r, gC, bC};
            fwrite(px, 1, 3, fp);
        }
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
        int s1 = seed_min + row;
        for (int col = 0; col < W; ++col)
        {
            int s2 = seed_min + col;

            mpz_t a, b;
            mpz_init_set_ui(a, (unsigned long)s1);
            mpz_init_set_ui(b, (unsigned long)s2);

            int excl = 0, terms = 0;

            /* term 0 */
            terms++;
            if ((s1 > 10) && (s1 % 2 == 0 || s1 % 5 == 0))
                excl++;

            /* term 1 */
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

            double frac = (terms > 0) ? ((double)excl / (double)terms) : 0.0; /* 0..1 */
            unsigned char r, g, bC;
            palette_fire(frac, &r, &g, &bC);
            unsigned char px[3] = {r, g, bC};
            fwrite(px, 1, 3, fp);
        }
    }
    fclose(fp);
    printf("[HEATMAP] Wrote %s (quick-exclusion fraction)\n", path);
}

/* ---------- Main ---------- */

int main(int argc, char **argv)
{
    /* CLI-compatible defaults (unchanged) */
    int seed_min = 1, seed_max = 50, window = 100, top = 25, max_terms = 200000, target_digits = 3800;
    const char *out_path = "results.txt";
    const char *batch_path = "batch-results.txt";

    for (int i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "--seed-min") && i + 1 < argc)
            seed_min = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--seed-max") && i + 1 < argc)
            seed_max = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--window") && i + 1 < argc)
            window = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--top") && i + 1 < argc)
            top = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--max-terms") && i + 1 < argc)
            max_terms = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--target-digits") && i + 1 < argc)
            target_digits = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--out") && i + 1 < argc)
            out_path = argv[++i];
        else if (!strcmp(argv[i], "--batch") && i + 1 < argc)
            batch_path = argv[++i];
    }

    printf("== FLO Prime Finder (C + GMP) ==\n");
    printf("seed_min=%d seed_max=%d window=%d top=%d target_digits=%d max_terms=%d\n",
           seed_min, seed_max, window, top, target_digits, max_terms);

    /* Stage 1: Seed scan with early-abandon diagnostics */
    int total = (seed_max - seed_min + 1) * (seed_max - seed_min + 1);
    SeedStat *stats = (SeedStat *)malloc(sizeof(SeedStat) * total);
    SeedStat *stats_grid = (SeedStat *)malloc(sizeof(SeedStat) * total); /* preserve pre-sort order for heatmap */
    int idx = 0;

    const double keepup_threshold = 0.80; /* abandon if cannot reach 80% of current best */
    int best_hits_so_far = 0;
    int abandoned_count = 0;

    clock_t t0 = clock();
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

            if (abandoned_flag)
                abandoned_count++;
            if (hits > best_hits_so_far)
                best_hits_so_far = hits;
        }
    }
    double elapsed = (double)(clock() - t0) / CLOCKS_PER_SEC;
    printf("[STATS] scanned %d seed pairs in %.2fs | early-abandoned=%d (threshold=%.2f of best)\n",
           total, elapsed, abandoned_count, keepup_threshold);

    /* Generate regular + illuminated heatmaps from the grid snapshot */
    write_heatmap_ppm(stats_grid, seed_min, seed_max, "heatmap.ppm");
    write_gcd_heatmap_ppm(seed_min, seed_max, "heatmap_gcd.ppm");
    write_exclusions_heatmap_ppm(seed_min, seed_max, window, "heatmap_exclusions.ppm");

    /* Additional diagnostics */
    int gcd_gt1 = 0, diagonal = 0;
    for (int i = 0; i < total; ++i)
    {
        int a = stats_grid[i].s1, b = stats_grid[i].s2, g;
        int x = a, y = b;
        while (y)
        {
            int t = x % y;
            x = y;
            y = t;
        }
        g = x;
        if (g > 1)
            gcd_gt1++;
        if (stats_grid[i].s1 == stats_grid[i].s2)
            diagonal++;
    }
    printf("[DIAG] gcd(s1,s2)>1 seeds: %d/%d | diagonal seeds: %d\n", gcd_gt1, total, diagonal);

    /* Leaderboard (unchanged) */
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
    top = 1;
    printf("[SELECT] Using seed with max prime hits (ties randomized): (%d,%d) hits=%d (tie-group=%d)\n",
           stats[0].s1, stats[0].s2, stats[0].hits, tie_count);

    /* Stage 2: Focused search (unchanged output) */
    FILE *out = fopen(out_path, "w");
    fprintf(out, "# Focused FLO prime search\n# target_digits=%d max_terms=%d\n", target_digits, max_terms);

    double expected_checks = expected_prime_trials(target_digits, /*odd_only=*/1, /*corrections=*/0);
    printf("Estimated primality tests per prime @ %d digits (odd-only): %.1f\n", target_digits, expected_checks);
    fprintf(out, "Estimated primality tests per prime @ %d digits (odd-only): %.1f\n", target_digits, expected_checks);

    mpz_t a, b;
    mpz_init(a);
    mpz_init(b);

    int total_found = 0;

    for (int i = 0; i < (top < total ? top : total); ++i)
    {
        int s1 = stats[i].s1, s2 = stats[i].s2;
        fprintf(out, "\n== Seed (%d,%d) rank=%d window_hits=%d density=%.3f\n",
                s1, s2, i + 1, stats[i].hits, stats[i].density);
        printf("[SEARCH] Seed (%d,%d)\n", s1, s2);

        mpz_set_ui(a, (unsigned long)s1);
        mpz_set_ui(b, (unsigned long)s2);

        /* warm-up to target digits */
        int terms_advanced = 0;
        while (mpz_sizeinbase(b, 10) < (unsigned)target_digits)
        {
            flo_next(a, b);
            terms_advanced++;
            if (terms_advanced > max_terms)
                break;
        }

        if (mpz_sizeinbase(b, 10) < (unsigned)target_digits)
        {
            printf("  Reached max warm-up terms without hitting %d digits. Skipping seed.\n", target_digits);
            fprintf(out, "  Reached max warm-up terms without hitting %d digits. Skipping.\n", target_digits);
            continue;
        }

        clock_t t_start = clock();
        long long checks = 0;
        int found_for_seed = 0;

        for (int t = 0; t < max_terms; ++t)
        {
            clock_t _chk_t0 = clock();

            if (mpz_divisible_ui_p(b, 2) || mpz_divisible_ui_p(b, 5))
            {
                flo_next(a, b);
                continue;
            }

            checks++;
            int pr = mpz_probab_prime_p(b, 25);

            clock_t _chk_t1 = clock();
            double _chk_sec = (double)(_chk_t1 - _chk_t0) / CLOCKS_PER_SEC;
            double _elapsed = (double)(_chk_t1 - t_start) / CLOCKS_PER_SEC;
            double _avg_chk = (checks > 0) ? (_elapsed / (double)checks) : _chk_sec;
            double _remain = (expected_checks > (double)checks) ? (expected_checks - (double)checks) : 0.0;
            double _eta_sec = _remain * _avg_chk;

            fprintf(stderr, "\rcheck #%lld | last=%.6fs avg=%.6fs | est. remain≈%.1f checks | ETA≈%.1fs   ",
                    checks, _chk_sec, _avg_chk, _remain, _eta_sec);

            if (pr > 0)
            {
                clock_t t_end = clock();
                double secs = (double)(t_end - t_start) / CLOCKS_PER_SEC;

                char *prime_str = mpz_get_str(NULL, 10, b);
                int dcount = (int)strlen(prime_str);

                double eff_ratio = (checks > 0) ? (expected_checks / (double)checks) : 0.0;

                printf("FOUND probable prime | seed=(%d,%d) digits=%d\n", s1, s2, dcount);
                printf("  time_to_find: %.3f sec | checks: %lld | expected_checks: %.1f | efficiency (expected/actual): %.3f\n",
                       secs, checks, expected_checks, eff_ratio);
                printf("  prime: %s\n", prime_str);

                fprintf(out, "FOUND probable prime | seed=(%d,%d) digits=%d\n", s1, s2, dcount);
                fprintf(out, "  time_to_find: %.3f sec | checks: %lld | expected_checks: %.1f | efficiency (expected/actual): %.3f\n",
                        secs, checks, expected_checks, eff_ratio);
                fprintf(out, "  prime: %s\n", prime_str);
                fflush(out);

                free(prime_str);
                total_found++;
                found_for_seed = 1;
                break;
            }

            flo_next(a, b);
        }

        if (!found_for_seed)
        {
            printf("  No prime found within %d checks-window (post-threshold) for this seed.\n", max_terms);
            fprintf(out, "  No prime found within %d checks-window (post-threshold).\n", max_terms);
        }
    }

    fprintf(out, "\n# total_found=%d\n", total_found);
    fclose(out);

    mpz_clear(a);
    mpz_clear(b);
    free(stats_grid);
    free(stats);
    return 0;
}
