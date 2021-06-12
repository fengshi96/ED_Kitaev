//
// Created by shifeng on 3/13/21.
//
#include "BLAS.h"
#include "Vector.h"

// ============================================================================
// =                           Vector Norm Wrapper                            =
// ============================================================================
double norm(const std::vector<dcomplex>& v){
    int n = v.size();
    int incx = 1;
    return BLAS::dznrm2_(&n, &(v[0]), &incx);
}

float norm(const std::vector<fcomplex>& v){
    int n = v.size();
    int incx = 1;
    return BLAS::scnrm2_(&n, &(v[0]), &incx);
}

double norm(const std::vector<double>& v){
    int n = v.size();
    int incx = 1;
    return BLAS::dnrm2_(&n, &(v[0]), &incx);
}

float norm(const std::vector<float>& v){
    int n = v.size();
    int incx = 1;
    return BLAS::snrm2_(&n, &(v[0]), &incx);
}

// ============================================================================
// =                           Vector Scalar Wrapper                          =
// ============================================================================
void vscal(const dcomplex& a, std::vector<dcomplex>& X) {
    int incx = 1;
    int n = X.size();
    BLAS::zscal_(&n, &a, &(X[0]), &incx);
}

void vscal(const fcomplex& a, std::vector<fcomplex>& X) {
    int incx = 1;
    int n = X.size();
    BLAS::cscal_(&n, &a, &(X[0]), &incx);
}

void vscal(const double& a, std::vector<double>& X) {
    int incx = 1;
    int n = X.size();
    BLAS::dscal_(&n, &a, &(X[0]), &incx);
}

void vscal(const float& a, std::vector<float>& X) {
    int incx = 1;
    int n = X.size();
    BLAS::sscal_(&n, &a, &(X[0]), &incx);
}

// ============================================================================
// =                Vector Vector Multiplication Wrapper                      =
// ============================================================================
// X dot Y
double dot(std::vector<double>& X, std::vector<double>& Y) {
    int incx = 1; int incy = 1;  // stride
    int n = X.size();
    return BLAS::ddot_(&n, &(X[0]), &incx, &(Y[0]), &incy);
}

float dot(std::vector<float>& X, std::vector<float>& Y) {
    int incx = 1; int incy = 1;  // stride
    int n = X.size();
    return BLAS::sdot_(&n, &(X[0]), &incx, &(Y[0]), &incy);
}

dcomplex dot(std::vector<dcomplex>& X, std::vector<dcomplex>& Y) {
    int incx = 1; int incy = 1;  // stride
    int n = X.size();
    return BLAS::zdotu_(&n, &(X[0]), &incx, &(Y[0]), &incy);
}

fcomplex dot(std::vector<fcomplex>& X, std::vector<fcomplex>& Y) {
    int incx = 1; int incy = 1;  // stride
    int n = X.size();
    return BLAS::cdotu_(&n, &(X[0]), &incx, &(Y[0]), &incy);
}

// X.conjugate dot Y
dcomplex cdot(std::vector<dcomplex>& X, std::vector<dcomplex>& Y) {
    int incx = 1; int incy = 1;  // stride
    int n = X.size();
    return BLAS::zdotc_(&n, &(X[0]), &incx, &(Y[0]), &incy);
}

fcomplex cdot(std::vector<fcomplex>& X, std::vector<fcomplex>& Y) {
    int incx = 1; int incy = 1;  // stride
    int n = X.size();
    return BLAS::cdotc_(&n, &(X[0]), &incx, &(Y[0]), &incy);
}
