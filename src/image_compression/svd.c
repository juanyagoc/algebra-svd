/* File: /src/image_compression/svd.c
 * Description: Main functions to obtain the SVD decomposition
 */

#include <string.h>
#include <math.h>
#include <stdio.h>

#include "qr_decomposition/tools.h"
#include "qr_decomposition/qr_decomp.h"

// TODO(#1): Improve the QR computations to it doest start to diverge when the number of iteriations on the
//       obtain_diagonal_matrix main loop is too high.

/**
 * This function constructs a given hermitian matrix M multiplying A by its transpose. M contains all 0 except possibly
 * an approximation of the singular values of A in its main diagonal. The process involves iterative steps in which
 * M collapses each time to a more diagonal form.
 *
 * @param A The given matrix of dimensions m by n
 * @param AT the transpose of A, n x m
 * @param M hermitian matrix, first A * A^t, then diagonal with the singuar values of A on it
 * @param m rows of A
 * @param n cols of A
 * @param iterations number of iterations, more implies better approximation of the s.v. of A.
 */
void obtain_diagonal_matrix(const double* A, double* AT, double* M, const int m, const int n, const int iterations)
{
    transpose(A, AT, m, n);
    matrix_mult(AT, A, M, n, n, m);

    for (int i = 0; i < iterations; i++) {
        double tau[n];
        double Q[n * n];
        double R[n * n];
        double RES[n * n];

        memset(Q, 0, sizeof(double) * n * n);
        memset(R, 0, sizeof(double) * n * n);
        memset(RES, 0, sizeof(double) * n * n);
        memset(tau, 0, sizeof(double) * n);

        build_qr_decomposition(M, tau, Q, R, n, n);

        matrix_mult(R, Q, RES, n, n, n);

        // Sum of the elements out of the diagonal, looking for the convergence time
        double off_diagonal_norm = 0.0;
        for (int k = 0; k < n; k++) {
            for (int j = 0; j < n; j++) {
                if (k != j) {
                    off_diagonal_norm += fabs(RES[k * n + j]);
                }
            }
        }

        memcpy(M, RES, sizeof(double) * n * n);

        if (off_diagonal_norm < 1e-3) { break; } // Stop when the convergence reaches a level of acceptance
    }
}

void obtain_right_singular_vectors(const double* A, double* AT, double* V, const int m, const int n, const int iterations)
{
    double M[n * n];
    memset(M, 0, sizeof(double) * n * n);
    transpose(A, AT, m, n);
    matrix_mult(AT, A, M, n, n, m);

    for (int i = 0; i < iterations; i++) {
        double tau[n];
        double Q[n * n];
        double R[n * n];
        memset(Q, 0, sizeof(double) * n * n);
        memset(R, 0, sizeof(double) * n * n);
        memset(tau, 0, sizeof(double) * n);

        build_qr_decomposition(M, tau, Q, R, n, n);
        matrix_mult(Q, R, M, n, n, n);

        if (i == iterations - 1) {
            memcpy(V, Q, sizeof(double) * n * n);
        }
    }
}

void obtain_left_singular_vectors(const double* A, const double* V, double* U, const double* singular_values, const int m, const int n)
{
    double inverse_sigma[n * n];
    memset(inverse_sigma, 0, sizeof(double) * n * n);

    for (int i = 0; i < n; i++) {
        if (singular_values[i] != 0) {
            inverse_sigma[i * n + i] = 1.0/singular_values[i];
        }
    }
    matrix_mult(V, inverse_sigma, inverse_sigma, n, n, n);
    matrix_mult(A, inverse_sigma, U, m, m, n);
}