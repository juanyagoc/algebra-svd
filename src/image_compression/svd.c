/* File: /src/image_compression/svd.c
 * Description: Main functions to obtain the SVD decomposition
 */

#include <string.h>

#include "qr_decomposition/tools.h"
#include "qr_decomposition/qr_decomp.h"

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
    for (int iter = 0; iter < iterations; iter++) {
        double tau[m];
        double Q[n * n];
        double R[n * n];
        double RES[n * n];
        memset(Q, 0, sizeof(double) * n * n);
        memset(R, 0, sizeof(double) * n * n);
        memset(RES, 0, sizeof(double) * n * n);
        memset(tau, 0, sizeof(double) * m);
        build_qr_decomposition(M, tau, Q, R, m, n);
        matrix_mult(Q, R, RES, m, n, n);
        for (int i = 0; i < m * n; i++) {
            M[i] = RES[i];
        }
    }
    print_matrix("Matriz Final: ", M, m, n);
}