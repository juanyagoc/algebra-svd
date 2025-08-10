/* File: /src/image_compression/main.c
 * Description: Main file.
 */

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "qr_decomposition/tools.h"
#include "svd.h"

void compile_svd(const double* A, double* U, double* M, double* V, const int m, const int n,
                 const int iterations) {
    double* AT = malloc(m * n * sizeof(double));
    memset(AT, 0, m * n * sizeof(double));
    obtain_diagonal_matrix(A, AT, M, m, n, iterations);
    obtain_right_singular_vectors(A, AT, V, m, n, iterations);
    const int minimum = min(m, n);
    double* singular_values = malloc(sizeof(double) * minimum);
    for (int i = 0; i < min(m, n); i++) {
        if (M) {
            const double M_i = M[i * n + i];
            singular_values[i] = sqrt(M_i);
        }
    }
    obtain_left_singular_vectors(A, V, U, singular_values, m, n);
}

int main() {
    const int m = 8; // Change the dimensions for testing other matrices
    const int n = 6;
    double* AT = malloc(m * n * sizeof(double));
    memset(AT, 0, m * n * sizeof(double));
    double* M = malloc(n * n * sizeof(double));
    double U[m * n];
    double V[m * n];

    const double A[8 * 6] = {
        7.0,  -2.0,   0.5,  3.0,  -1.0,  4.0,   2.0,  -3.0,
       -1.0,   5.0,   2.0, -4.0,   6.0, -2.0,   1.0,   3.0,
        3.0,   0.0,  -2.0,  5.0,  -1.0,  7.0,  -3.0,   2.0,
        0.0,  -6.0,   4.0,  1.0,   2.0, -5.0,   3.0,  -2.0,
        8.0,  -3.0,   1.0, -2.0,   0.0,  4.0,  -1.0,   6.0,
       -4.0,   2.0,   5.0, -1.0,   3.0,  0.0,  -2.0,   7.0
    };

    const double B[6] = {
        3,2,
        2,3,
        2,-2
    };

    compile_svd(A, U, M, V, m,  n, 20);
    print_matrix("Matriz diagoal M", M, n, n);
    print_matrix("Matriz V", V, n, n);
    print_matrix("Matriz U", U, m, m);
    return 0;
}
