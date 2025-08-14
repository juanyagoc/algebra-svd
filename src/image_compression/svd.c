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

        const double a = M[(n-2)*n + (n-2)];
        const double b = M[(n-1)*n + (n-2)];
        const double c = M[(n-1)*n + (n-1)];
        const double delta = (a - c) / 2.0;
        const double sign_delta = (delta >= 0.0) ? 1.0 : -1.0;
        const double mu = c - (sign_delta * b * b) / (fabs(delta) + sqrt(delta * delta + b * b));

        for (int k = 0; k < n; k++) {
            M[k*n + k] -= mu;
        }


        build_qr_decomposition(M, tau, Q, R, n, n);

        matrix_mult(R, Q, RES, n, n, n);

        for (int k = 0; k < n; k++) {
            RES[k * n + k] += mu;
        }

        memcpy(M, RES, sizeof(double) * n * n);

        // Sum of the elements out of the diagonal, looking for the convergence time
        double off_diagonal_norm = 0.0;
        for (int k = 0; k < n; k++) {
            for (int j = 0; j < n; j++) {
                if (k != j) {
                    off_diagonal_norm += fabs(RES[k * n + j]);
                }
            }
        }

        if (off_diagonal_norm < 1e-15) { break; } // Stop when the convergence reaches a level
                                                  // of acceptance
    }
}


/**
 * This method computes the right matrix V of the SVD decomposition using the QR decomposition
 * iteratively applied to M = A^T * A as in the previous method, in fact they coud be combined
 * aboiding repetitive calculations this is not done in my implementation.
 * On each iteration the next matrix M is equal to R * Q and this time we store the rotations
 * that represent Q on the matrix V, so V_{k} = V_{k-1} * Q.
 *
 * @param A Input matrix of size m x n remember we work in column-major.
 * @param AT Buffer for keeping the transpose of A it is important to assure its preallocated
 *           when the method is called.
 * @param V Output matrix, will be size n x n. Starts as the identity and ends being the orthogonal
 *          approximation of the singuar vectors of A.
 * @param m Number of rows of A
 * @param n number of columns of A
 * @param iterations number of maximum iterations to perform on the call normally reducing the
 *                   approximation error.
 *
 * @complexity O(m * n^2 + n^3 + iterations)
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
 * With this method we can obtain the last matrix in the SVD decomposition. The general notation
 * used is A = U \Sigma V, so knowing the original matrix A, the diagonal matrix \Sigma and V
 * the only thing lwft is to compute U with A V^-1 \Sigma^-1. Because \Sigma is a diagonal matrix
 * its inverse consists of 1/{the entries on the diagonal} and regarding V we know its orthogonal
 * so its inverse is its transpose, multiplying those tree we obtain U.
 *
 * @param A               Input matrix m x n
 * @param V               Orthogonal matrix of right singular vectors (n x n)
 * @param U               Output matrix m x n containing the left singular vectors
 * @param singular_values Array containing the singular values
 * @param m               Number of rows of A
 * @param n               Number of columns of A.
 */

void obtain_left_singular_vectors(const double* A, const double* V, double* U,
                                  const double* singular_values, const int m, const int n)
{
    double inverse_sigma[n * n];
    memset(inverse_sigma, 0, sizeof(double) * n * n);

    for (int i = 0; i < n; i++) {
        if (singular_values[i] > 1e-15) {
            inverse_sigma[i * n + i] = 1.0 / singular_values[i];
        }
    }

    double inverse_v[n * n];
    transpose(V, inverse_v, n, n);

    double temp[n * n];
    matrix_mult(inverse_v, inverse_sigma, temp, n, n, n);
    matrix_mult(A, temp, U, m, n, n);
}