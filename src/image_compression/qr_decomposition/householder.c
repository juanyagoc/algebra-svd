/* File: src/image_compression/qr_decomposition/householder.c
 * Description: This file computes the reflection householder vector and its corresponding matrix
 */

#include "householder.h"
#include "tools.h"
#include <stdio.h>
#include <math.h>

/**
 * TODO(#1): Implement the blocking QR decomposition version as seen in the Golub book.
 */

/**
 * Householder reflecion computation based on the algorithm described in the Golub book.
 * It receives a vector 'x' and reflects it on other computed vector so the reflection 'v' lands on
 * a scalar multiple of a orthonormal vector of the form {1, 0, ..., 0} with the desired dimension.
 * This reflectoin can be used to construct the Householder matrix H wich will be applied to A.
 *
 * @param x the given vector to be transformed by a reflection
 * @param v the output vector being the Householder vector
 * @param tau the scalar factor used in the matrix computation
 * @param n the dimension of the vector 'x'
 */
void householder_vector(const double* x, double* v, double* tau, const int n)
{
    if (n <= 0) { printf("Error en la longitud n = %d", n); return; }
    double sigma = norm(x + 1, n - 1);
    sigma = sigma * sigma; // Precission limit 1e+-160 risk of underflow/overflow

    if (sigma == 0.0) {
        *tau = 0.0;
        if (x[0] >= 0) {
            *tau = 0.0;
        } else {
            *tau = -2.0;
        }
        v[0] = 1.0;
        for (int i = 1; i < n; i++) {
            v[i] = 0.0;
        }
        return;
    }

    for (int i = 1; i < n; i++) {
        v[i] = x[i];
    }

    const double mu = sqrt(x[0] * x[0] + sigma);

    if (x[0] <= 0.0) {
        v[0] = x[0] - mu;
    } else {
        v[0] = -sigma / (x[0] + mu);
    }

    const double v0_sqr = v[0] * v[0];
    *tau = 2.0 * v0_sqr / (sigma + v0_sqr);

    for (int i = 1; i < n; i++) {
        v[i] = v[i] / v[0];
    }
    v[0] = 1.0;
}

/**
 * @brief Applies the Householder reflections to the given matrix A.
 *
 * This method iteratively transforms the matrix A to obtaing the upper triangular matrix R
 * of the QR decomposition. Although A is not exactly R because the lower triangular halt of A stores
 * the necessary vectors to form the matrix Q.
 *
 * @param A Pointer to the given matrix A with m rows and n columns stored in column-major order as used in LAPACK
 *          library, when the method ends A results as described below.
 * @param tau Array of length n which stores the scalars needed for each Householder vector
 * @param m number of rows of A
 * @param n number of columns of A
 */
void compute_householder_matrices(double* A, double* tau, const int m, const int n)
{
    for (int i = 0; i < n && i < m - 1; i++) {
        int const len = m - i;
        double v[len];
        householder_vector(&A[i + i * m], v, &tau[i], len);

        for (int j = i; j < n; j++) {
            double dot = 0.0;
            for (int k = 0; k < len; k++) {
                dot += v[k] * A[i + k + j * m];
            }
            dot *= tau[i];
            for (int k = 0; k < len; k++) {
                A[i + k + j * m] -= dot * v[k];
            }
        }

        for (int k = 1; k < len; k++) {
            A[i + k + i * m] = v[k];
        }
    }
}
