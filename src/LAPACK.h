//
// Created by shifeng on 3/11/21.
//

#ifndef CPPKIT_LAPACK_H
#define CPPKIT_LAPACK_H
#include <complex>

typedef std::complex<double> dcomplex;
typedef std::complex<float> fcomplex;

namespace LAPACK {
    extern "C" {
        // ============================================================================
        // =                             Diagonalization                              =
        // ============================================================================
        void zheev_(char *, char *, int *, dcomplex *, int *, double *, dcomplex *, int *, double *, int *);
        void cheev_(char *, char *, int *, fcomplex *, int *, float *, fcomplex *, int *, float *, int *);
        void dsyev_(char *, char *, int *, double *, int *, double *, double *, int *, int *);
        void ssyev_(char *, char *, int *, float *, int *, float *, float *, int *, int *);
        void dgeev_(char *, char *, int *, double *, int *, double *, double *, double *, int *, double *, int *,
                    double *, int *, int *);
        void sgeev_(char *, char *, int *, float *, int *, float *, float *, float *, int *, float *, int *, float *,
                    int *, int *);
        void zgeev_(char *, char *, int *, dcomplex *, int *, dcomplex *, dcomplex *, int *, dcomplex *,
                    int *, dcomplex *, int *, double *, int *);
        void cgeev_(char *, char *, int *, fcomplex *, int *, fcomplex *, fcomplex *, int *, fcomplex *,
                    int *, fcomplex *, int *, float *, int *);
    }
}

#endif //CPPKIT_LAPACK_H
