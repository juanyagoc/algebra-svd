/* File: src/image_compression/qr_decomposition/qr_decomp.c
 * Description: This file contains the functions that form the resulting Q and R matrices
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "householder.h"

/**
 * Performs the operations as seen in the algorithm 5.1.5 of Golubs book, obtaining each Q_k from the previous ones
 * by the operation Q_k = Q_{k-1} - tau * v * v^t * Q_{k-1}. So it basically takes the dot product of v and Q,
 * multiplyes it by tau and v, then substracts that from Q, assigning that new value to Q.
 *
 * With this implementation we avoid creating the transformation matrix I - tau * v * v^t by applying the operations
 * directly on Q.
 *
 * @param Q Unitary square matrix dimenions m x m, store in column-major as all this algorithms always do.
 * @param v Householder vector of length m - i. Its initial value is always 1.0 for convention.
 * @param tau corresponding scalar coefficient for the v vector
 * @param m rows of A, also rows and columns of Q
 * @param i the dimensions of the active/current submatrix where we are applying the transformation.
 */
void apply_householder_transform(double* Q, const double* v, const double tau, const int m, const int i)
{
    const int length = m - i;

    for (int j = 0; j < m; j++) {

        double *q_ptr = Q + i + j * m;
        double dot = 0.0;
        const double *vk = v;
        double *qk = q_ptr;

        for (int k = 0; k < length; ++k) {
            dot += (*vk++) * (*qk++);
        }

        dot *= tau;

        if (dot != 0.0) {
            vk = v;
            qk = q_ptr;
            for (int k = 0; k < length; ++k) {
                qk[k] -= dot * vk[k];
            }
        }
    }
}

/**
 * @brief Given a matrix A obtains Q and R such that A = QR.
 *
 * This functions stracts the R matrix from the upper triangular portion of A after A its being transformed by
 * 'compute_householder_matrices' function. Then inicializes Q as the identity matrix. Then takes each column lower
 * part below the main diagonal and stores it on v which has always the first value 1.0, then calls the
 * 'apply_householder_tarnsform' function with the apropiated tau value corresponting to the current v to build Q.
 *
 * @param A the matrix after being transformed, dimensions m by n storen in column-major order.
 * @param tau Array of doubles that stores the scalars associated with each Householder reflection
 * @param Q Output the resulting unitary square matrix of dimensions m x m
 * @param R Output the resulting upper triangular matrix of dimensions m x n
 * @param m number of rows of A and Q
 * @param n number of columns of A and R
 */
void build_qr_decomposition(double* A, double* tau, double* Q, double* R, const int m, const int n)
{
    compute_householder_matrices(A, tau, m, n);

    memset(R, 0, sizeof(double) * m * n);
    for (int i = 0; i < n; ++i) {
        const int max_row = (i < m) ? i : (m - 1); // cap the row
        for (int j = 0; j <= max_row; ++j) {
            R[j + i * m] = A[j + i * m];
        }
    }

    memset(Q, 0, sizeof(double) * m * m);
    for (int i = 0; i < m; i++) {
        Q[i + i * m] = 1.0;
    }

    double *v = malloc(sizeof(double) * m);
    if (!v) { return; }

    for (int i = n - 1; i >= 0; i--) {
        const int length = m - i;
        if (length <= 0) { continue; }
        const double tau_i = tau[i];

        if (fabs(tau_i) < 1e-15) { continue; }

        v[0] = 1.0;
        for (int j = 1; j < length; j++) {
            v[j] = A[(i + j) + i * m];
        }

        apply_householder_transform(Q, v, tau_i, m, i);
    }

    free(v);
}