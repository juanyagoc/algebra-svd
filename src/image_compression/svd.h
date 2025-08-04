#pragma once
#ifdef __cplusplus
extern "C" {
#endif

    void obtain_diagonal_matrix(const double* A, double* AT, double* M, int m, int n, int iterations);

    void obtain_right_singular_vectors(const double* A, double* AT, double* V, int m, int n, int iterations);

    void obtain_left_singular_vectors(const double* A, const double* V, double* U, const double* singular_values, int m, int n);

#ifdef __cplusplus
}
#endif