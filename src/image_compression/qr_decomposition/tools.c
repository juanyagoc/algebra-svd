#include <stdio.h>
#include "tools.h"
#include <math.h>

// always using row-major matrices
// 2-norm with double pass technique for possible overflow
double norm(const double* x, const int n) { // n is the lenght of vector x
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

double inner_product(const double* x, const double* y, const int n) {
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

// Better if implements 'blocking' for big matrices, also implement matrix-vector multiplication function
void matrix_mult(const double* A, const double* B, double* C, int m, int n, int l) {
    // C[m x x] = A[m x l] * B[l * n]
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

void print_matrix(const char* name, const double* mat, int rows, int cols) {
    printf("%s =\n", name);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            printf("%8.3f ", mat[i + j * cols]);  // column-major
        }
        printf("\n");
    }
    printf("\n");
}
