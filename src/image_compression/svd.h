#pragma once
#ifdef __cplusplus
extern "C" {
#endif

void obtain_diagonal_matrix(const double* A, double* AT, double* M, int m, int n, int iterations);

#ifdef __cplusplus
}
#endif