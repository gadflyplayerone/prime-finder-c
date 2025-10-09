/*
 * Prime Finder via Fibonacci‑Like Operator (C + GMP)
 *
 * Build: see Makefile
 */

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef struct {
    int s1, s2;
    int hits;
    double density;
} SeedStat;

static int is_small_prime(unsigned long n) {
    static const int small[] = {2,3,5,7,11,13,17,19,23,29,31,37};
    size_t i;
    if (n < 2) return 0;
    for (i=0;i<sizeof(small)/sizeof(small[0]);++i) {
        if (n == (unsigned long)small[i]) return 1;
        if (n % small[i] == 0) return 0;
    }
    return 1; // maybe prime
}

static void flo_next(mpz_t a, mpz_t b) {
    // (a,b) -> (b, a+b)
    mpz_t tmp;
    mpz_init(tmp);
    mpz_add(tmp, a, b);
    mpz_set(a, b);
    mpz_set(b, tmp);
    mpz_clear(tmp);
}

static int count_primes_in_window(int s1, int s2, int window) {
    mpz_t a,b;
    mpz_init_set_ui(a, (unsigned long)s1);
    mpz_init_set_ui(b, (unsigned long)s2);
    int hits = 0;

    // term 0
    if (s1 > 10 && (s1%2==0 || s1%5==0)) { /* skip */ }
    else {
        if (mpz_probab_prime_p(a, 10) > 0) hits++;
    }

    // term 1
    if (s2 > 10 && (s2%2==0 || s2%5==0)) { /* skip */ }
    else {
        if (mpz_probab_prime_p(b, 10) > 0) hits++;
    }

    for (int i=2;i<window;i++) {
        flo_next(a,b);
        if (mpz_sizeinbase(b, 2) > 4 && mpz_divisible_ui_p(b, 2)) continue;
        if (mpz_divisible_ui_p(b, 5)) continue;
        if (mpz_probab_prime_p(b, 10) > 0) hits++;
    }
    mpz_clear(a); mpz_clear(b);
    return hits;
}

static int cmp_seedstat(const void* A, const void* B) {
    const SeedStat* a = (const SeedStat*)A;
    const SeedStat* b = (const SeedStat*)B;
    if (a->hits != b->hits) return (b->hits - a->hits);
    if (a->density < b->density) return 1;
    if (a->density > b->density) return -1;
    return 0;
}

int main(int argc, char** argv) {
    int seed_min=1, seed_max=50, window=100, top=25, max_terms=200000, target_digits=3800;
    const char* out_path = "results.txt";
    const char* batch_path = "batch-results.txt";

    for (int i=1;i<argc;i++) {
        if (!strcmp(argv[i], "--seed-min") && i+1<argc) seed_min = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--seed-max") && i+1<argc) seed_max = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--window") && i+1<argc) window = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--top") && i+1<argc) top = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--max-terms") && i+1<argc) max_terms = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--target-digits") && i+1<argc) target_digits = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--out") && i+1<argc) out_path = argv[++i];
        else if (!strcmp(argv[i], "--batch") && i+1<argc) batch_path = argv[++i];
    }

    printf("== FLO Prime Finder (C + GMP) ==\n");
    printf("seed_min=%d seed_max=%d window=%d top=%d target_digits=%d max_terms=%d\n",
           seed_min, seed_max, window, top, target_digits, max_terms);

    int total = (seed_max - seed_min + 1) * (seed_max - seed_min + 1);
    SeedStat* stats = (SeedStat*)malloc(sizeof(SeedStat)*total);
    int idx = 0;

    clock_t t0 = clock();
    for (int s1=seed_min; s1<=seed_max; ++s1) {
        for (int s2=seed_min; s2<=seed_max; ++s2) {
            int hits = count_primes_in_window(s1, s2, window);
            stats[idx].s1 = s1;
            stats[idx].s2 = s2;
            stats[idx].hits = hits;
            stats[idx].density = (double)hits / (double)(window>0?window:1);
            idx++;
        }
    }
    double elapsed = (double)(clock()-t0)/CLOCKS_PER_SEC;
    printf("[STATS] scanned %d seed pairs in %.2fs\n", total, elapsed);

    qsort(stats, total, sizeof(SeedStat), cmp_seedstat);

    FILE* batch = fopen(batch_path, "w");
    fprintf(batch, "# Seed Scan Leaderboard\n# window=%d\n", window);
    fprintf(batch, "rank,s1,s2,hits,density\n");
    for (int i=0;i< (top*2<total?top*2:total); ++i) {
        fprintf(batch, "%d,%d,%d,%d,%.3f\n", i+1, stats[i].s1, stats[i].s2, stats[i].hits, stats[i].density);
        printf("%d,%d,%d,%d,%.3f\n", i+1, stats[i].s1, stats[i].s2, stats[i].hits, stats[i].density);
    }
    fclose(batch);

    // Focused search
    FILE* out = fopen(out_path, "w");
    fprintf(out, "# Focused FLO prime search\n# target_digits=%d max_terms=%d\n", target_digits, max_terms);

    mpz_t a,b;
    mpz_init(a); mpz_init(b);

    int found = 0;
    for (int i=0;i< (top<total?top:total); ++i) {
        int s1 = stats[i].s1, s2 = stats[i].s2;
        fprintf(out, "\n== Seed (%d,%d) rank=%d window_hits=%d density=%.3f\n",
                s1, s2, i+1, stats[i].hits, stats[i].density);
        printf("[SEARCH] Seed (%d,%d) …\n", s1, s2);

        mpz_set_ui(a, (unsigned long)s1);
        mpz_set_ui(b, (unsigned long)s2);

        // Iterate
        for (int term=0; term<max_terms; ++term) {
            // ensure 'b' is current term
            if (term>=2) {
                mpz_add(b, a, b); // b = a + b
                mpz_swap(a, b);   // (a,b) = (b,a) then next add will produce next term in 'b'
                mpz_add(b, a, b); // produce current b again due to swap logic
                mpz_swap(a, b);
                // The above double-step keeps (a,b) aligned; simpler alternative is a tmp func (omitted for speed).
            }
            int bits = mpz_sizeinbase(b, 10);
            if (bits < target_digits) continue;
            if (mpz_divisible_ui_p(b, 2) || mpz_divisible_ui_p(b, 5)) continue;
            int pr = mpz_probab_prime_p(b, 25);
            if (pr > 0) {
                found++;
                char* s = mpz_get_str(NULL, 10, b);
                int dcount = (int)strlen(s);
                fprintf(out, "FOUND probable prime | seed=(%d,%d) term=%d digits=%d\n", s1, s2, term, dcount);
                printf("FOUND probable prime | seed=(%d,%d) term=%d digits=%d\n", s1, s2, term, dcount);
                free(s);
                fflush(out);
            }
        }
        printf("[DONE] Seed (%d,%d)\n", s1, s2);
    }

    fprintf(out, "\n# total_found=%d\n", found);
    fclose(out);
    mpz_clear(a); mpz_clear(b);
    free(stats);
    return 0;
}