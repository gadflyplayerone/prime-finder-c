#include "flo_predict.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

static const double LOG10_PHI = 0.20898764024997873;

static const struct { int m, r; double lift; } DEFAULT_ARMS[] = {
    {24,11,2.182540}, {21,2,2.025463}, {12,11,1.981707}, {30,23,1.909722},
    {27,11,1.827485}, {20,3,1.805556}, {18,5,1.760913}, {28,11,1.736111}
};

typedef struct { int i; double s; } Pair;

static int cmp_pair_desc(const void *A, const void *B) {
    const Pair *a = (const Pair*)A;
    const Pair *b = (const Pair*)B;
    if (a->s < b->s) return 1;
    if (a->s > b->s) return -1;
    return 0;
}

int flo_estimate_index_for_digits(int digits) {
    return (int)ceil(digits / LOG10_PHI);
}

void flo_predictor_init(FloPredictor *P, int target_digits, int window_half_width) {
    memset(P, 0, sizeof(*P));
    size_t N = sizeof(DEFAULT_ARMS)/sizeof(DEFAULT_ARMS[0]);
    P->arms = (ResidueArm*)calloc(N, sizeof(ResidueArm));
    P->n_arms = N;
    for (size_t i=0;i<N;i++) {
        P->arms[i].m = DEFAULT_ARMS[i].m;
        P->arms[i].r = DEFAULT_ARMS[i].r;
        P->arms[i].lift = DEFAULT_ARMS[i].lift;
        P->arms[i].alpha = 1.0;
        P->arms[i].beta  = 1.0;
    }
    P->target_digits = target_digits;
    P->n0_estimate = flo_estimate_index_for_digits(target_digits);
    P->window_start = P->n0_estimate - window_half_width;
    P->window_end   = P->n0_estimate + window_half_width;

    srand((unsigned int)time(NULL));
}

void flo_predictor_destroy(FloPredictor *P) {
    if (P->arms) free(P->arms);
    memset(P, 0, sizeof(*P));
}

static double bandit_mean(double a, double b) {
    return a / (a + b);
}

double flo_index_score(const FloPredictor *P, int index) {
    double s = 0.0;
    for (size_t k=0;k<P->n_arms;k++) {
        const ResidueArm *a = &P->arms[k];
        if (index % a->m == a->r) {
            s += a->lift + bandit_mean(a->alpha, a->beta);
        }
    }
    return s;
}

size_t flo_rank_window(const FloPredictor *P, int start, int end, int *out_indices, size_t max_out) {
    if (end < start) return 0;
    int len = end - start + 1;
    Pair *buf = (Pair*)malloc(sizeof(Pair) * (size_t)len);
    for (int i=0;i<len;i++) {
        buf[i].i = start + i;
        buf[i].s = flo_index_score(P, start + i);
    }
    qsort(buf, (size_t)len, sizeof(Pair), cmp_pair_desc);
    int n = (len < (int)max_out) ? len : (int)max_out;
    for (int k=0;k<n;k++) out_indices[k] = buf[k].i;
    free(buf);
    return (size_t)n;
}

void flo_bandit_update(FloPredictor *P, int index, int is_prime) {
    for (size_t k=0;k<P->n_arms;k++) {
        ResidueArm *a = &P->arms[k];
        if (index % a->m == a->r) {
            if (is_prime) a->alpha += 1.0;
            else          a->beta  += 1.0;
        }
    }
}
