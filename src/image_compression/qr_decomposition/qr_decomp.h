#pragma once
#ifdef __cplusplus
extend C {}
#endif

    void apply_householder_transform(double* Q, const double* v, double tau, int m, int i);

    void build_qr_decomposition(const double* A, const double* tau, double* Q, double* R, int m, int n);

#ifdef __cplusplus
}
#endif