#include <stdlib.h>
#include "householder.h"
#include "tools.h"
#include <stdio.h>

int main() {
    const int m = 3, n = 3;
    double tau[m];

    double A[9] = {
        12, 6, -4,
       -51, 167, 24,
         4, -68, -41
    };

    compute_householder_matrices(A, tau, m, n);

    print_matrix("A", A, m, n);

    return 0;
}