#include "tools.h"
#include "householder.h"
#include "qr_decomp.h"

int main() {
    const int m = 3;
    const int n = 3;
    const int block_size = 1;
    double tau[m];

    // Remember to define 'A' transposed because we are working column-major
    double A[9] = {
        12, 6, -4,
       -51, 167, 24,
         4, -68, -41
    };

    // Noral 'A' without transposing it
    double B[9] = {
        12, -51, 4,
        6, 167, -68,
        -4, 24, -41
    };

    double Q[m * m], R[m * n];
    compute_householder_matrices(A, tau, m, n);
    build_qr_decomposition(A, tau, Q, R, m, n);

    print_matrix("Q", Q, m, m);
    print_matrix("R", R, m, n);

    return 0;
}