#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "qr_decomposition/tools.h"
#include "svd.h"

int main() {
    const int m = 2;
    const int n = 3;
    double tau[m];
    double* AT = malloc(m * n * sizeof(double));
    memset(AT, 0, m * n * sizeof(double));
    double* M = malloc(n * n * sizeof(double));
    double U[m * n];
    double V[m * n];
    double singular_values[2]= {0};

    // Remember to define 'A' transposed because we are working column-major
    const double A[9] = {
        12, 6, -4,
       -51, 167, 24,
         4, -68, -41
    };

    const double B[6] = {
        3,2,
        2,3,
        2,-2
    };

    obtain_diagonal_matrix(B, AT, M, m, n, 15);
    print_matrix("M", M, n, n);

    for (int i = 0; i < min(m, n); i++) {
        singular_values[i] = sqrt(M[i * n + i]);
    }

    //print_matrix("Valores singulares: ", singular_values, 1, min(m, n));

    //obtain_right_singular_vectors(B, AT, V, m, n, 15);
    //print_matrix("Matriz V: ", V, n, n);
    //obtain_left_singular_vectors(B, V, U, singular_values, m, n);
    //print_matrix("Matriz U: ", U, m, m);

    return 0;
}



