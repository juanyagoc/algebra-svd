#include "../../../src/image_compression/qr_decomposition/tools.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define TOLERANCE 1e-10

int comparation_with_tolerance(const double a, const double b) {
    return fabs(a - b) < TOLERANCE;
}

void test_norm() {
    const double x[3] = {1.0, 0.5, 6.0};
    const double real_norm = sqrt(37.25);
    const double norm_x = norm(x, 3);

    if (!comparation_with_tolerance(real_norm, norm_x)) {
        printf("ERROR: The 'norm' output is not the correct. Spected %f, got %f\n",
                real_norm, norm_x);
        exit(1);
    }

    const double zero_vector[5] = {0};
    if (norm(zero_vector, 5) != 0) {
        printf("ERROR: The 'norm' function with zero vector doest output 0.0\n");
        exit(1);
    }

    const double scalar = 3.45;
    const double scaled_x[3] = { x[0] * scalar, x[1] * scalar, x[2] * scalar };

    const double norm_scaled = norm(scaled_x, 3);
    if (!comparation_with_tolerance(norm_scaled, fabs(scalar) * norm_x)) {
        printf("ERROR: Scalar property failed. Expected %f, got %f\n",
                fabs(scalar) * norm_x, norm_scaled);
        exit(1);
    }

    printf("§ test_norm OK\n");
}

int main(void) {
    test_norm();
    printf("All tests passed\n");
    return 0;
}