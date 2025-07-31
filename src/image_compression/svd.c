/* File: /src/image_compression/svd.c
 * Description: Main functions to obtain the SVD decomposition
 */

#include "qr_decomposition/tools.h"

void obtain_hermitian_matrix(const double* A, double* AT, double* M, const int m, const int n)
{
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            AT[j + i * n] = A[i + j * m];
        }
    }
    matrix_mult(A, AT, M, n, n, m);
}