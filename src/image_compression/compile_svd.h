#pragma once
#ifdef __cplusplus
extern "C" {
#endif

    void compile_svd(const double* A, double* U, double* M, double* V, const int m, const int n,
                     int iterations);

#ifdef __cplusplus
}
#endif