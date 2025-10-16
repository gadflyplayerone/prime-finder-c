#ifndef FLO_PREDICT_H
#define FLO_PREDICT_H

#include <stddef.h>

typedef struct {
    int m;
    int r;
    double lift;
    double alpha;
    double beta;
} ResidueArm;

typedef struct {
    ResidueArm *arms;
    size_t n_arms;
    int window_start;
    int window_end;
    int target_digits;
    int n0_estimate;
} FloPredictor;

void flo_predictor_init(FloPredictor *P, int target_digits, int window_half_width);
void flo_predictor_destroy(FloPredictor *P);
int  flo_estimate_index_for_digits(int digits);
double flo_index_score(const FloPredictor *P, int index);
size_t flo_rank_window(const FloPredictor *P, int start, int end, int *out_indices, size_t max_out);
void flo_bandit_update(FloPredictor *P, int index, int is_prime);

#endif
