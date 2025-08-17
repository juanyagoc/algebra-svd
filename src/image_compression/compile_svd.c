/* File: /src/image_compression/compile_svd.c
 * Description: Here the SVD matrices get constructed.
 */

#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "svd.h"

/**
 *
 *
 * @param A
 * @param U
 * @param M
 * @param V
 * @param m
 * @param n
 * @param iterations
 */
void compile_svd(const double* A, double* U, double* M, double* V, const int m, const int n,
                 const int iterations) {
    double* AT = malloc(m * n * sizeof(double));
    memset(AT, 0, m * n * sizeof(double));
    obtain_diagonal_matrix(A, AT, M, m, n, iterations);
    obtain_right_singular_vectors(A, AT, V, m, n, iterations);

    double* singular_values = malloc(sizeof(double) * n);

    for (int i = 0; i < n; i++) {
        const double d = M[i * n + i];
        singular_values[i] = sqrt(d > 0.0 ? d : 0.0);
    }

    int* idx = malloc(sizeof(int) * n);
    for (int i = 0; i < n; i++) idx[i] = i;
    for (int i = 0; i < n - 1; i++) {
        int best = i;
        for (int j = i + 1; j < n; j++) {
            if (singular_values[j] > singular_values[best]) best = j;
        }
        if (best != i) {
            const double tmpv = singular_values[i];
            singular_values[i] = singular_values[best];
            singular_values[best] = tmpv;
            const int tmpi = idx[i];
            idx[i] = idx[best];
            idx[best] = tmpi;
        }
    }

    double* M_sorted = malloc(sizeof(double) * n * n);
    memset(M_sorted, 0, sizeof(double) * n * n);
    for (int k = 0; k < n; k++) {
        const double sigma_k = singular_values[k];
        M_sorted[k * n + k] = sigma_k * sigma_k;
    }
    memcpy(M, M_sorted, sizeof(double) * n * n);

    obtain_left_singular_vectors(A, V, U, singular_values, m, n);

    free(M_sorted);
    free(idx);
    free(singular_values);
    free(AT);
}
