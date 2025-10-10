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

static int count_primes_in_window(int s1, int s2, int window)
{
    mpz_t a, b;
    mpz_init_set_ui(a, (unsigned long)s1);
    mpz_init_set_ui(b, (unsigned long)s2);
    int hits = 0;

    // term 0
    if (!(s1 > 10 && (s1 % 2 == 0 || s1 % 5 == 0)))
    {
        if (mpz_probab_prime_p(a, 10) > 0)
            hits++;
    }
    // term 1
    if (!(s2 > 10 && (s2 % 2 == 0 || s2 % 5 == 0)))
    {
        if (mpz_probab_prime_p(b, 10) > 0)
            hits++;
    }
    for (int i = 2; i < window; i++)
    {
        flo_next(a, b);
        if (mpz_divisible_ui_p(b, 2) || mpz_divisible_ui_p(b, 5))
            continue;
        if (mpz_probab_prime_p(b, 10) > 0)
            hits++;
    }
    mpz_clear(a);
    mpz_clear(b);
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

    /* Stage 1: Seed scan (unchanged) */
    int total = (seed_max - seed_min + 1) * (seed_max - seed_min + 1);
    SeedStat *stats = (SeedStat *)malloc(sizeof(SeedStat) * total);
    int idx = 0;

    clock_t t0 = clock();
    for (int s1 = seed_min; s1 <= seed_max; ++s1)
    {
        for (int s2 = seed_min; s2 <= seed_max; ++s2)
        {
            int hits = count_primes_in_window(s1, s2, window);
            stats[idx].s1 = s1;
            stats[idx].s2 = s2;
            stats[idx].hits = hits;
            stats[idx].density = (double)hits / (double)(window > 0 ? window : 1);
            idx++;
        }
    }
    double elapsed = (double)(clock() - t0) / CLOCKS_PER_SEC;
    printf("[STATS] scanned %d seed pairs in %.2fs\n", total, elapsed);

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

    /* ---- NEW: choose seed pair with MAX hits; randomize among ties ---- */
    srand((unsigned)time(NULL));
    int max_hits = stats[0].hits; /* after sort, stats[0] has the maximum hits */
    int tie_count = 0;
    while (tie_count < total && stats[tie_count].hits == max_hits)
    {
        tie_count++;
    }
    int chosen_idx = (tie_count > 0) ? (rand() % tie_count) : 0;

    /* Move chosen seed into stats[0] to reuse downstream logic with minimal changes */
    if (chosen_idx != 0)
    {
        SeedStat tmp = stats[0];
        stats[0] = stats[chosen_idx];
        stats[chosen_idx] = tmp;
    }
    /* Force single-seed focused search */
    top = 1;
    printf("[SELECT] Using seed with max prime hits (ties randomized): (%d,%d) hits=%d (tie-group=%d)\n",
           stats[0].s1, stats[0].s2, stats[0].hits, tie_count);

    /* Stage 2: Focused search with the five requested behaviors */
    FILE *out = fopen(out_path, "w");
    fprintf(out, "# Focused FLO prime search\n# target_digits=%d max_terms=%d\n", target_digits, max_terms);

    /* Pre-compute expected checks for reporting (odd-only assumption) */
    double expected_checks = expected_prime_trials(target_digits, /*odd_only=*/1, /*corrections=*/0);
    printf("Estimated primality tests per prime @ %d digits (odd-only): %.1f\n",
           target_digits, expected_checks);
    fprintf(out, "Estimated primality tests per prime @ %d digits (odd-only): %.1f\n",
            target_digits, expected_checks);

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

        /* Initialize sequence */
        mpz_set_ui(a, (unsigned long)s1);
        mpz_set_ui(b, (unsigned long)s2);

        /* 1) Warm-up: advance without any prime checks until target digit length is reached */
        int terms_advanced = 0;
        while (mpz_sizeinbase(b, 10) < (unsigned)target_digits)
        {
            flo_next(a, b);
            terms_advanced++;
            if (terms_advanced > max_terms)
                break; /* safety */
        }

        if (mpz_sizeinbase(b, 10) < (unsigned)target_digits)
        {
            printf("  Reached max warm-up terms without hitting %d digits. Skipping seed.\n", target_digits);
            fprintf(out, "  Reached max warm-up terms without hitting %d digits. Skipping.\n", target_digits);
            continue;
        }

        /* 2–3) Start timing and count only checks AFTER threshold is met */
        clock_t t_start = clock();
        long long checks = 0;
        int found_for_seed = 0;

        for (int t = 0; t < max_terms; ++t)
        {
            /* --- per-check timing start --- */
            clock_t _chk_t0 = clock();

            /* Quick filters do NOT count as checks */
            if (mpz_divisible_ui_p(b, 2) || mpz_divisible_ui_p(b, 5))
            {
                flo_next(a, b);
                continue;
            }

            /* 3) Count this as a primality check */
            checks++;
            int pr = mpz_probab_prime_p(b, 25);

            /* --- per-check timing end + ETA --- */
            clock_t _chk_t1 = clock();
            double _chk_sec = (double)(_chk_t1 - _chk_t0) / CLOCKS_PER_SEC;          /* duration of THIS prime_check */
            double _elapsed = (double)(_chk_t1 - t_start) / CLOCKS_PER_SEC;          /* total timed runtime since threshold */
            double _avg_chk = (checks > 0) ? (_elapsed / (double)checks) : _chk_sec; /* rolling avg per check */
            double _remain = (expected_checks > (double)checks) ?                    /* est. checks remaining */
                                 (expected_checks - (double)checks)
                                                                : 0.0;
            double _eta_sec = _remain * _avg_chk; /* ETA to first prime (remaining) */

            fprintf(stderr,
                    "\rcheck #%lld | last=%.6fs avg=%.6fs | est. remain≈%.1f checks | ETA≈%.1fs   ",
                    checks, _chk_sec, _avg_chk, _remain, _eta_sec);

            if (pr > 0)
            {
                /* 2) Stop timer on first hit */
                clock_t t_end = clock();
                double secs = (double)(t_end - t_start) / CLOCKS_PER_SEC;

                char *prime_str = mpz_get_str(NULL, 10, b);
                int dcount = (int)strlen(prime_str);

                /* 4) Report stats vs expected + 5) Print prime */
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
                break; /* stop on first prime for this seed */
            }

            /* advance sequence */
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
    free(stats);
    return 0;
}
