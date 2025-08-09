/* File: /src/image_compression/svd.c
 * Description: Main functions to obtain the SVD decomposition
 */

#include <string.h>
#include <math.h>
#include <stdio.h>

#include "qr_decomposition/tools.h"
#include "qr_decomposition/qr_decomp.h"

/*
 * TODO(#1): Improve the QR computations to it doest start to diverge when the number of
 *           iteriations on the obtain_diagonal_matrix main loop is too high. Better convergence
 *           with 'Wilkinson shif'
 */

/**
 * This function constructs a given hermitian matrix M multiplying A by its transpose.
 * M contains all 0 except possibly an approximation of the singular values of A in its
 * main diagonal. The process involves iterative steps in which M collapses each time to
 * a more diagonal form.
 *
 * @param A The given matrix of dimensions m by n
 * @param AT the transpose of A, n x m
 * @param M hermitian matrix, first A * A^t, then diagonal with the singuar values of A on it
 * @param m rows of A
 * @param n cols of A
 * @param iterations number of iterations, more implies better approximation of the s.v. of A.
 */
void obtain_diagonal_matrix(const double* A, double* AT, double* M, const int m, const int n,
                            const int iterations)
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

        if (off_diagonal_norm < 1e-3) { break; } // Stop when the convergence reaches a level
                                                 // of acceptance
    }
}

/**
 *
 * @param A
 * @param AT
 * @param V
 * @param m
 * @param n
 * @param iterations
 */
void obtain_right_singular_vectors(const double* A, double* AT, double* V, const int m,
                                   const int n, const int iterations)
{
    double tau[n];
    double Q[n * n];
    double R[n * n];
    double M[n * n];
    double PREV_M[n * n];
    double NEXT_V[n * n];

    memset(M, 0, sizeof(double) * n * n);
    transpose(A, AT, m, n); // AT está bien calculada
    matrix_mult(AT, A, M, n, n, m); // M está bien calculada, el error no es aquí

    memset(V, 0, sizeof(double) * n * n);
    for (int i = 0; i < n; i++) {
        V[i + i * n] = 1.0;
    }

    for (int i = 0; i < iterations; i++) {
        memset(Q, 0, sizeof(double) * n * n);
        memset(R, 0, sizeof(double) * n * n);
        memset(PREV_M, 0, sizeof(double) * n * n);
        memset(NEXT_V, 0, sizeof(double) * n * n);
        memset(tau, 0, sizeof(double) * n);

        build_qr_decomposition(M, tau, Q, R, n, n);
        //print_matrix("TAU", tau, n, 1);
        matrix_mult(R, Q, PREV_M, n, n, n);
        memcpy(M, PREV_M, sizeof(double) * n * n);

        matrix_mult(V, Q, NEXT_V, n, n, n);
        memcpy(V, NEXT_V, sizeof(double) * n * n);
    }
}

/**
 *
 * @param A
 * @param V
 * @param U
 * @param singular_values
 * @param m
 * @param n
 */
void obtain_left_singular_vectors(const double* A, const double* V, double* U,
                                  const double* singular_values, const int m, const int n)
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