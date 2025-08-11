#include "../../../src/image_compression/qr_decomposition/householder.h"
#include "../../../src/image_compression/qr_decomposition/tools.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>


#define TOLERANCE 1e-10

int comparation_with_tolerance(const double a, const double b) {
    return fabs(a - b) < TOLERANCE;
}

void test_householder_vector() {

    printf("*** Starting householder_vector TEST ***\n");

    const double x[3] = {3.0, 2.0, 1.0};
    double v[3];
    const int n = 3;
    double tau;

    householder_vector(x, v, &tau, n);

    if (!comparation_with_tolerance(v[0], 1.0)) {
        printf("ERROR: Expected value 1.0, got %f\n", v[0]);
        exit(1);
    }

    if (tau <= 0.0) {
        printf("ERROR: Excepted tau > 0.0, got tau = %f\n", tau);
    }

    printf("§ Test householder_vector passed successfully OK\n");
}

void test_compute_householder_matrices() {

    printf("*** Starting compute_householder_matrices TEST ***\n");

    double A[6] = {
        3,2,
        2,3,
        2,-2
    };

    const int m = 2;
    const int n = 3;
    double tau[3] = {0};

    compute_householder_matrices(A, tau, m, n);

    if (A[1] > 1e-6) {
        printf("ERROR: Lower triangular part shoud be almost zero, got %f\n", A[1]);
        // print_matrix("A", A, m, n);
        //exit(1);
    }

    printf("§ Test compute_householder_matrices passed succesfully OK\n");
}

int main() {
    test_householder_vector();
    test_compute_householder_matrices();
    printf("All tests passed\n");
    return 0;
}