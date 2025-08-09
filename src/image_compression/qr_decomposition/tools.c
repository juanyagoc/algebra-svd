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
    if (n <= 0) { return 0.0; }
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
 * Implements the inner product operation between two given vector x and y which have to have the same length n.
 * Also preventing overflow with the double pass operations.
 *
 * @param x Fisrt vector
 * @param y secont vector
 * @param n length of both vectors.
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
 * Given a matrix A and a possibly empty matrix AT, transposes A and saves the result into AT.
 * Note that if A is given as matrix to transpose and as matrix to save the transposition it does notperforms
 * the operation in place the reuslt woud be a matrix A with the lower triangular part copied from its upper
 * triangular part.
 *
 * @param A Input matrix to stranspose
 * @param AT transpose saving argument
 * @param m number of rows of A, columns of AT
 * @param n number fo columns of AT, rows of A.
 */
void transpose(const double* A, double * AT, const int m, const int n)
{
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            AT[j + i * n] = A[i + j * m];
        }
    }
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
 * @param l number of columns of A and rows of B, necessary to be the same.
 */
void matrix_mult(const double* A, const double* B, double* C, const int m, const int n, const int l) {
    for (int i = 0; i < m * n; i++) {
        C[i] = 0.0;
    }

    // Column-major order multiplication
    for (int j = 0; j < n; ++j) {
        for (int k = 0; k < l; ++k) {
            for (int i = 0; i < m; ++i) {
                C[j * m + i] += A[k * m + i] * B[j * l + k];
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
 * @param cols number of columns of the matrix.
 */
void print_matrix(const char* name, const double* mat, const int rows, const int cols)
{
    printf("%s =\n", name);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            printf("%8.4f ", mat[j * rows + i]);  // column-major
        }
        printf("\n");
    }
    printf("\n");
}

int min(const int a, const int b) {
    return (a < b) ? a : b;
}
