/* File: /src/image_compression/svd.c
 * Description: Main functions to obtain the SVD decomposition
 */

#include "qr_decomposition/tools.h"
#include "qr_decomposition/qr_decomp.h"

void obtain_diagonal_matrix(const double* A, double* AT, double* M, const int m, const int n)
{
    transpose(A, AT, m, n);
    matrix_mult(A, AT, M, n, n, m);
    for (int iter = 0; iter < 10; iter++) {
        double tau[m];
        double Q[n * n];
        double R[n * n];
        double M_copy[n * n];

        for (int i = 0; i < n * n; ++i) {
            M_copy[i] = M[i];
        }

        build_qr_decomposition(M_copy, tau, Q, R, m, n);
        matrix_mult(R, Q, M, m, n, n);
    }

    print_matrix("M", M, m, n);
}