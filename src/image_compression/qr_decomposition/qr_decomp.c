#include <stdlib.h>
#include <string.h>
#include "tools.h"

/**
 * Applies the householder trasnformactions.
 *
 * @param Q unitary matrix
 * @param v vectors
 * @param tau scalars
 * @param m dimensions
 * @param i number
 */
void apply_householder_transform(double* Q, const double* v, const double tau, const int m, const int i) {
    for (int j = 0; j < m; j++) {
        double dot_product = 0.0;
        for (int k = 0; k < m - i; ++k) {
            dot_product += v[k] * Q[i + k + j * m];
        }
        dot_product *= tau;
        for (int k = 0; k < m - i; k++) {
            Q[i + k + j * m] -= dot_product * v[k];
        }
    }
}

/**
 *
 *
 * @param A
 * @param tau
 * @param Q
 * @param R
 * @param m
 * @param n
 */
void build_qr_decomposition(const double* A, const double* tau, double* Q, double* R, const int m, const int n) {
    memset(R, 0, m * n * sizeof(double));
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            R[i + j * m] = A[i + j * m];
        }
    }

    memset(Q, 0, m * m * sizeof(double));
    for (int i = 0; i < m; i++) {
        Q[i + i * m] = 1.0;
    }

    for (int i = n - 1; i >= 0; i--) {
        double v[m];
        memset(v, 0, m * sizeof(double));
        v[0] = 1.0;
        const int len = m - i;
        for (int j = 1; j < len; j++) {
            v[j] = A[i + j + i * m];
        }
        apply_householder_transform(Q, v, tau[i], m, i);
    }
}