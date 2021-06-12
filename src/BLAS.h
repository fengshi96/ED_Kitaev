//
// Created by shifeng on 3/11/21.
//

#ifndef CPPKIT_BLAS_H
#define CPPKIT_BLAS_H

#include <complex>
typedef std::complex<double> dcomplex;
typedef std::complex<float> fcomplex;

namespace BLAS {
    extern "C" {
        // ============================================================================
        // =    Level 1 BLAS          DOT: dot product of vectors                     =
        // ============================================================================
        // X dot Y
        double ddot_(const int*, const double *, const int* ,const double*, const int* );
        float sdot_(const int*, const float *, const int*, const float*, const int* );
        dcomplex zdotu_(const int*, const dcomplex *, const int*, const dcomplex*, const int* );
        fcomplex cdotu_(const int*, const fcomplex *, const int*, const fcomplex*, const int* );

        // X.conj dot Y
        dcomplex zdotc_(const int*, const dcomplex *, const int*, const dcomplex*, const int* );
        fcomplex cdotc_(const int*, const fcomplex *, const int*, const fcomplex*, const int* );

        // ****************************************************************************
        // *                          SCAL: rescale vectors
        // ****************************************************************************
        void sscal_(const int *, const float *, float *, const int *);
        void dscal_(const int *, const double *, double *, const int *);
        void cscal_(const int *, const fcomplex *, fcomplex *, const int *);
        void zscal_(const int *, const dcomplex *, dcomplex *, const int *);

        // ****************************************************************************
        // *                          NRM2: vector norm
        // ****************************************************************************
        float scnrm2_(const int *, const fcomplex *, const int *);
        double dznrm2_(const int *, const dcomplex *, const int *);
        double dnrm2_(const int *, const double *, const int *);
        float snrm2_(const int *, const float *, const int *);


        // ============================================================================
        // =    Level 2 BLAS          GEMV: matrix vector multiply                    =
        // ============================================================================
        void sgemv_(char*,int*,int*,const float*,const float*,int*,
                               const float*,int*,const float*,float*,int*);
        void dgemv_(char*,int*,int*,const double*,const double*,int*,
                               const double*,int*,const double*,double*,int*);
        void cgemv_(char*,int*,int*,const fcomplex*,const fcomplex*,
                               int*,const fcomplex*,int*,const fcomplex*, fcomplex*,int*);
        void zgemv_(char*,int*,int*,const dcomplex*,const dcomplex*,
                               int*,const dcomplex*,int*,const dcomplex*,dcomplex*,int*);
        // ****************************************************************************
        // *                          HEMV: hermitian matrix vector multiply
        // ****************************************************************************
        void chemv_(char*,int*,const fcomplex*,const fcomplex*,int*,const fcomplex*,
                                int*,const fcomplex*,fcomplex*,int*);
        void zhemv_(char*,int*,const dcomplex*,const dcomplex*,
                                int*,const dcomplex*,int*,const dcomplex*,dcomplex*,int*);


        // ============================================================================
        // = Level 3 BLAS             GEMM: matrix matrix multiply                    =
        // ============================================================================
        void sgemm_(char*,char*,int*,int*,int*,const float*,
                               const float*,int*,const float*,int*,const float*,float*,int*);
        void dgemm_(char*,char*,int*,int*,int*,const double*,
                               const double*,int*,const double*,int*,const double*,double*,int*);
        void cgemm_(char*,char*,int*,int*,int*,const fcomplex*,
                               const fcomplex*,int*,const fcomplex*,int*,const fcomplex*,fcomplex*,int*);
        void zgemm_(char*,char*,int*,int*,int*,const dcomplex*,
                               const dcomplex*,int*,const dcomplex*,int*,const dcomplex*,dcomplex*,int*);
    }
}

#endif //CPPKIT_BLAS_H
