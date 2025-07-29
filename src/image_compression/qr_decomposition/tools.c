/* File: src/image_compression/qr_decomposition/tools.c
 * Description: This file contains basic matrix and vector operations
 */

#include <stdio.h>
#include <math.h>

#include "tools.h"

/**
 * Computes the 2-norm of a given vector. Implements the double pass technique preventing possible overflow,
 * always working in column-major order.
 *
 * @param x the vector
 * @param n the length of the vector x
 * @return double representing the norm of the vector.
 */
double norm(const double* x, const int n)
{
    double max = 0.0;
    for (int i = 0; i < n; i++) {
        const double abs_x = fabs(x[i]);
        if (max < abs_x) {
            max = abs_x;
        }
    }

    if (max == 0) { return 0.0; }

    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        const double sca_x = x[i] / max;
        sum += sca_x * sca_x;
    }
    return max * sqrt(sum);
}

/**
 *
 * @param x
 * @param y
 * @param n
 * @return
 */
double inner_product(const double* x, const double* y, const int n)
{
    double max = 0.0;
    for (int i = 0; i < n; i++) {
        const double abs_x = fabs(x[i]);
        const double abs_y = fabs(y[i]);
        if (max < abs_x) {
            max = abs_x;
        }
        if (max < abs_y) {
            max = abs_y;
        }
    }
    if (max == 0) { return 0.0; }
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        const double sca_x = x[i] / max;
        const double sca_y = y[i] / max;
        sum += sca_x * sca_y;
    }
    return sum;
}

/**
 * Multiplies two given matrices A, B and saves the result on the matrix C. The orden of the iterations in the 'for'
 * loops is the recomended in the Golub and BanLoan book for better memory access.
 *
 * @param A First matrix given in column-major order with dimensions m x l
 * @param B second matrix, with dimensions l x n
 * @param C resulting matrix with dimensions m x n
 * @param m number of row of A
 * @param n number of columns of B
 * @param l number of columns of A and rows of B, necessary to be the same
 */
void matrix_mult(const double* A, const double* B, double* C, const int m, const int n, const int l)
{
    // bucle order recomended in 'Matrix Computations' by Golub for better memory running
    for (int i = 0; i < m; ++i) {
        for (int k = 0; k < l; ++k) {
            const double a_ik = A[i * l + k];
            for (int j = 0; j < n; ++j) {
                C[i * n + j] += a_ik * B[k * n + j];
            }
        }
    }
}

/**
 * Auxiliar method for printing matrices. Use it to test if the functions are working correctly by printing the
 * results on the terminal. Note it iterates through the matix in column-major order, it is necessary because all
 * the functions of this program work and store matrices in this order.
 *
 * @param name The name of the matrix it is going to print
 * @param mat the matrix itself
 * @param rows number of rows of the matrix
 * @param cols number of columns of the matrix
 */
void print_matrix(const char* name, const double* mat, const int rows, const int cols)
{
    printf("%s =\n", name);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            printf("%8.4f ", mat[i + j * cols]);  // column-major
        }
        printf("\n");
    }
    printf("\n");
}
