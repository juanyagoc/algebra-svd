#include "qr_decomposition/tools.h"
#include "qr_decomposition/qr_decomp.h"
#include "svd.h"

int main() {
    const int m = 3;
    const int n = 3;
    double tau[m];
    double AT[9];
    double M[9];

    // Remember to define 'A' transposed because we are working column-major
    const double A[9] = {
        12, 6, -4,
       -51, 167, 24,
         4, -68, -41
    };

    obtain_diagonal_matrix(A, AT, M, m, n);

    return 0;
}



