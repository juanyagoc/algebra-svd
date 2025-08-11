#include "../../../src/image_compression/qr_decomposition/qr_decomp.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define TOLERANCE 1e-10

int comparation_with_tolerance(const double a, const double b) {
    return fabs(a - b) < TOLERANCE;
}

static void set_identity(double* Q, const int m) {
    for (int j = 0; j < m; ++j) {
        for (int i = 0; i < m; ++i) {
            Q[i + j * m] = (i == j) ? 1.0 : 0.0;
        }
    }
}

static void build_embedded_householder(double* Hfull, const double* v, const double tau,
                                       const int m, const int i) {
    const int len = m - i;

    set_identity(Hfull, m);

    if (len <= 0) return;

    for (int r = 0; r < len; ++r) {
        for (int c = 0; c < len; ++c) {
            const double val = (r == c ? 1.0 : 0.0) - tau * v[r] * v[c];
            Hfull[(i + r) + (i + c) * m] = val;
        }
    }
}

static void test_apply_once(const int m, const int i, const double* v, const double tau) {
    double* Q = (double*)malloc(sizeof(double) * m * m);
    double* Hfull = (double*)malloc(sizeof(double) * m * m);
    if (!Q || !Hfull) {
        fprintf(stderr, "ERROR: memoria insuficiente\n");
        exit(1);
    }

    set_identity(Q, m);
    build_embedded_householder(Hfull, v, tau, m, i);

    apply_householder_transform(Q, v, tau, m, i);

    for (int col = 0; col < m; ++col) {
        for (int row = 0; row < m; ++row) {
            const double got = Q[row + col * m];
            const double exp = Hfull[row + col * m];
            if (!comparation_with_tolerance(got, exp)) {
                fprintf(stderr,
                        "ERROR: mismatch en (%d,%d): got=%+.12f, expected=%+.12f\n",
                        row, col, got, exp);
                free(Q);
                free(Hfull);
                exit(1);
            }
        }
    }

    free(Q);
    free(Hfull);
}

static void test_case_full_matrix() {
    printf("*** Starting full_matrix TEST ***\n");

    const int m = 4;
    const int i = 0;
    const double v[4] = {1.0, 0.2, -0.3, 0.5};
    double vtv = 0.0;
    for (int k = 0; k < m - i; ++k) vtv += v[k] * v[k];
    const double tau = 2.0 / vtv;

    test_apply_once(m, i, v, tau);
    printf("§ test_case_full_matrix OK\n");
}

static void test_case_submatrix() {
    printf("*** Starting case_submatrix TEST ***\n");

    const int m = 5;
    const int i = 2;
    const double v[3] = {1.0, -0.25, 0.4};
    double vtv = 0.0;
    for (int k = 0; k < m - i; ++k) {
        vtv += v[k] * v[k];
    }
    const double tau = 2.0 / vtv;

    test_apply_once(m, i, v, tau);
    printf("§ test_case_submatrix OK\n");
}

static void test_case_tau_zero() {
    printf("*** Starting case_tau_zero TEST ***\n");

    const int m = 3;
    const int i = 1;
    const double v[2] = {1.0, 0.7};
    const double tau = 0.0;
    double* Q = (double*)malloc(sizeof(double) * m * m);
    if (!Q) {
        fprintf(stderr, "ERROR: memory error creating Q\n");
        exit(1);
    }

    set_identity(Q, m);
    apply_householder_transform(Q, v, tau, m, i);

    for (int col = 0; col < m; ++col) {
        for (int row = 0; row < m; ++row) {
            const double exp = (row == col) ? 1.0 : 0.0;
            if (!comparation_with_tolerance(Q[row + col * m], exp)) {
                fprintf(stderr,
                        "ERROR: con tau=0, la matriz debe quedar identidad. Mismatch en (%d,%d)\n",
                        row, col );
                free(Q);
                exit(1);
            }
        }
    }

    free(Q);
    printf("§ test_case_tau_zero OK\n");
}

int main(void) {
    test_case_full_matrix();
    test_case_submatrix();
    test_case_tau_zero();
    printf("All tests passed\n");
    return 0;
}
