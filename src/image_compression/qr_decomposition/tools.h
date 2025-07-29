#pragma once
#ifdef __cplusplus
extern "C" {
#endif

    double norm(const double* x, int n);

    double inner_product(const double* x, const double* y, int n);

    void matrix_mult(const double* A, const double* B, double* C, int m, int n, int l);

    void print_matrix(const char* name, const double* mat, int rows, int cols);

#ifdef __cplusplus
}
#endif