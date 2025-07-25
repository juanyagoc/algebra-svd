#include <string.h>
#include <stdlib.h>

// Multiplica dos matrices cuadradas de tamaño m x m: C = A * B
void mat_mult(const double* A, const double* B, double* C, int m) {
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < m; ++j) {
            C[i * m + j] = 0.0;
            for (int k = 0; k < m; ++k)
                C[i * m + j] += A[i * m + k] * B[k * m + j];
        }
}

void build_QR_from_householder(const double* H, const double* A, double* Q, double* R, const int m, const int n, const int num_householders) {
    // Inicializa Q como identidad
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < m; ++j)
            Q[i * m + j] = (i == j) ? 1.0 : 0.0;

    // Acumula las transformaciones de Householder en Q (Q = H_0 * H_1 * ... * H_{k-1})
    double* temp = malloc(m * m * sizeof(double));
    for (int k = 0; k < num_householders; ++k) {
        const double* H_k = &H[k * m * m];
        mat_mult(H_k, Q, temp, m);
        memcpy(Q, temp, m * m * sizeof(double));
    }
    free(temp);

    // Calcula R = Q^T * A
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            R[i * n + j] = 0.0;
            for (int k = 0; k < m; ++k)
                R[i * n + j] += Q[k * m + i] * A[k * n + j];
        }
}