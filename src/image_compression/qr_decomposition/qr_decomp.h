#ifndef QR_DECOMP_H
#define QR_DECOMP_H

void build_QR_from_householder(const double* H, const double* A, double* Q, double* R, int m, int n, int num_householders);

#endif