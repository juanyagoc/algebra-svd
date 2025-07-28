#pragma once
#ifdef __cplusplus
extern "C" {
#endif

    // Computes the reflexion on given vector
    void householder_vector(const double* x, double* v, double* tau, int n);

    // Computes the Householder matrices H_k for a given input matrix A
    void compute_householder_matrices(double* A, double* tau, const int m, const int n);

    // Same computation but implementing the blocking method as describen in Golub's book
    void compute_blocked_qr(double* A, double* tau, const int m, const int n, const int b);

#ifdef __cplusplus
}
#endif