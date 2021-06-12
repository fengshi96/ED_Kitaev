#include "Matrix.h"
#include "LAPACK.h"
#include "BLAS.h"

// ============================================================================
// =                       Diagonalization Wrapper                            =
// ============================================================================
// diag double precision Hermitian matrix m
void diag(Matrix<dcomplex>& m, Vector<double>& evals, char option){
    char jobz=option;  // 'N':  Compute eigenvalues only; 'V':  Compute eigenvalues and eigenvectors.
    char uplo='U';  // 'U':  Upper triangle of A is stored;
    int n=m.rows();
    int lda=m.cols();
#ifdef DEBUG
    assert(m.IsHermitian());
#endif
    evals.resize(n);
    Vector<dcomplex> work(3);
    Vector<double> rwork(3*n-2);
    int info,lwork= -1;  // If LWORK = -1, then a workspace query is assumed

    evals.fill(0);

    // query:
    LAPACK::zheev_(&jobz,&uplo,&n,&(m(0,0)),&lda,&(evals[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    if (info!=0) {
        std::cerr<<"info="<<info<<" during query \n";
        perror("diag: zheev: failed with info!=0.\n");
    }
    lwork = static_cast<int>(std::real(work[0]))+1;
    work.resize(lwork+1);

    // real work:
    LAPACK::zheev_(&jobz,&uplo,&n,&(m(0,0)),&lda,&(evals[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    if (info!=0) {
        std::cerr<<"work info="<<info<<"\n";
        perror("diag: zheev: failed with info!=0.\n");
    }
    //for(int i=0;i<n;i++) { cout << eval[i] << " \t ";} cout << endl;
    //sort(eigs_.begin(),eigs_.end()); // sort Eigenvalues and Hamiltonian
}

// diag single precision Hermitian matrix m
void diag(Matrix<std::complex<float>> &m, std::vector<float>& evals, char option)
{
    char jobz=option;
    char uplo='U';
    int n=m.rows();
    int lda=m.cols();
#ifdef DEBUG
    assert(m.IsHermitian());
#endif
    std::vector<std::complex<float> > work(3);
    std::vector<float> rwork(3*n);
    int info,lwork= -1;

    evals.resize(n);

    // query:
    LAPACK::cheev_(&jobz,&uplo,&n,&(m(0,0)),&lda,&(evals[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    if (info!=0) {
        std::cerr<<"info="<<info<<" during query\n";
        throw std::runtime_error("diag: cheev_: failed with info!=0.\n");
    }

    const int NB = 256;
    lwork = std::max(1 + static_cast<int>(std::real(work[0])), (NB + 2)*n);
    work.resize(lwork);

    // real work:
    LAPACK::cheev_(&jobz,&uplo,&n,&(m(0,0)),&lda,&(evals[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    if (info!=0) {
        std::cerr<<"work info="<<info<<"\n";
        throw std::runtime_error("diag: cheev: failed with info!=0.\n");
    }
}


// diag double precision real symmetric matrix
void diag(Matrix<double> &m, std::vector<double>& evals, char option)
{
    char jobz=option;
    char uplo='U';
    int n=m.rows();
    int lda=m.cols();
#ifdef DEBUG
    assert(m.IsSymmetric());
#endif
    std::vector<double> work(3);
    int info;
    int lwork= -1;

    if (lda<=0) throw std::runtime_error("lda<=0\n");

    evals.resize(n);

    // query:
    LAPACK::dsyev_(&jobz,&uplo,&n,&(m(0,0)),&lda, &(evals[0]),&(work[0]),&lwork, &info);
    if (info!=0) {
        std::cerr<<"info="<<info<<" during query\n";
        throw std::runtime_error("diag: dsyev_: failed with info!=0.\n");
    }

    const int NB = 256;
    lwork = std::max(1 + static_cast<int>(work[0]), (NB + 2)*n);
    work.resize(lwork);
    // real work:
    LAPACK::dsyev_(&jobz,&uplo,&n,&(m(0,0)),&lda, &(evals[0]),&(work[0]),&lwork, &info);
    if (info!=0) {
        std::cerr<<"work info="<<info<<"\n";
        throw std::runtime_error("diag: dsyev_: failed with info!=0.\n");
    }
}

// diag single precision real symmetric matrix
void diag(Matrix<float> &m, std::vector<float>& evals,char option)
{
    char jobz=option;
    char uplo='U';
    int n=m.rows();
    int lda=m.cols();
#ifdef DEBUG
    assert(m.IsSymmetric());
#endif
    std::vector<float> work(3);
    int info;
    int lwork= -1;

    if (lda<=0) throw std::runtime_error("lda<=0\n");

    evals.resize(n);

    // query:
    LAPACK::ssyev_(&jobz,&uplo,&n,&(m(0,0)),&lda, &(evals[0]),&(work[0]),&lwork, &info);
    if (info!=0) {
        std::cerr<<"info="<<info<<" during query\n";
        throw std::runtime_error("diag: ssyev_: failed with info!=0.\n");
    }

    const int NB = 256;
    lwork = std::max(1 + static_cast<int>(work[0]), (NB + 2)*n);
    work.resize(lwork);

    // real work:
    LAPACK::ssyev_(&jobz,&uplo,&n,&(m(0,0)),&lda, &(evals[0]),&(work[0]),&lwork, &info);
    if (info!=0) {
        std::cerr<<"work info="<<info<<"\n";
        throw std::runtime_error("diag: ssyev_: failed with info!=0.\n");
    }
}

// diag double precision real non-symmetric matrix
void diag(Matrix<double> &m, std::vector<double>& evalsRe, std::vector<double>& evalsIm,
          std::vector<double> vr, char option){
    char jobvl = 'N';
    char jobvr = option;
    int n=m.rows();
    int lda=m.cols();
    assert(m.IsSquare());

    std::vector<double> vl;
    vr.resize(n * n);

    int ldvl = 1;
    int ldvr = n;
    std::vector<double> work(3);
    int lwork= -1;
    int info;

    evalsRe.resize(n);
    evalsIm.resize(n);

    // query:
    LAPACK::dgeev_(&jobvl, &jobvr, &n, &(m(0,0)), &lda, &(evalsRe[0]), &(evalsIm[0]),
           &(vl[0]), &ldvl, &(vr[0]), &ldvr, &(work[0]), &lwork, &info);
    if (info!=0) {
        std::cerr<<"info="<<info<<" during query\n";
        throw std::runtime_error("diag: dgeev_: failed with info!=0.\n");
    }
    lwork = static_cast<int>(work[0])+1;
    work.resize(lwork);

    // real work:
    LAPACK::dgeev_(&jobvl, &jobvr, &n, &(m(0,0)), &lda, &(evalsRe[0]), &(evalsIm[0]),
           &(vl[0]), &ldvl, &(vr[0]), &ldvr, &(work[0]), &lwork, &info);
    if (info!=0) {
        std::cerr<<"work info="<<info<<"\n";
        throw std::runtime_error("diag: dgeev_: failed with info!=0.\n");
    }
}

// diag single precision real non-symmetric matrix
void diag(Matrix<float> &m, std::vector<float>& evalsRe, std::vector<float>& evalsIm,
          std::vector<float> vr, char option){
    char jobvl = 'N';
    char jobvr = option;
    int n=m.rows();
    int lda=m.cols();
    assert(m.IsSquare());

    std::vector<float> vl;
    vr.resize(n * n);

    int ldvl = 1;
    int ldvr = n;
    std::vector<float> work(3);
    int lwork= -1;
    int info;

    evalsRe.resize(n);
    evalsIm.resize(n);

    // query:
    LAPACK::sgeev_(&jobvl, &jobvr, &n, &(m(0,0)), &lda, &(evalsRe[0]), &(evalsIm[0]),
           &(vl[0]), &ldvl, &(vr[0]), &ldvr, &(work[0]), &lwork, &info);
    if (info!=0) {
        std::cerr<<"info="<<info<<" during query\n";
        throw std::runtime_error("diag: sgeev_: failed with info!=0.\n");
    }
    lwork = static_cast<int>(work[0])+1;
    work.resize(lwork);

    // real work:
    LAPACK::sgeev_(&jobvl, &jobvr, &n, &(m(0,0)), &lda, &(evalsRe[0]), &(evalsIm[0]),
           &(vl[0]), &ldvl, &(vr[0]), &ldvr, &(work[0]), &lwork, &info);
    if (info!=0) {
        std::cerr<<"work info="<<info<<"\n";
        throw std::runtime_error("diag: sgeev_: failed with info!=0.\n");
    }
}

// diag double precision complex non-symmetric matrix
void diag(Matrix<dcomplex> &m, std::vector<dcomplex>& evals, std::vector<dcomplex>& vr, char option){
    char jobvl = 'N';
    char jobvr = option;
    int n=m.rows();
    int lda=m.cols();
    assert(m.IsSquare());

    int ldvl = 1;
    int ldvr = n;
    std::vector<dcomplex> vl;
    evals.resize(n);
    vr.resize(ldvr * n);
    std::vector<dcomplex> work(3);
    std::vector<double> rwork(2 * n);
    int lwork= -1;
    int info;

    // query:
    LAPACK::zgeev_(&jobvl, &jobvr, &n, &(m(0,0)), &lda, &(evals[0]), &(vl[0]), &ldvl, &(vr[0]), &ldvr,
           &(work[0]), &lwork, &(rwork[0]), &info);
    if (info!=0) {
        std::cerr<<"info="<<info<<" during query\n";
        throw std::runtime_error("diag: zgeev_: failed with info!=0.\n");
    }
    const int NB = 256;
    lwork = std::max(1 + static_cast<int>(std::real(work[0])), (NB + 2)*n);
    work.resize(lwork);

    // real work:
    LAPACK::zgeev_(&jobvl, &jobvr, &n, &(m(0,0)), &lda, &(evals[0]), &(vl[0]), &ldvl, &(vr[0]), &ldvr,
           &(work[0]), &lwork, &(rwork[0]), &info);
    if (info!=0) {
        std::cerr<<"work info="<<info<<"\n";
        throw std::runtime_error("diag: zgeev_: failed with info!=0.\n");
    }
}


// diag single precision complex non-symmetric matrix
void diag(Matrix<fcomplex> &m, std::vector<fcomplex>& evals, std::vector<fcomplex>& vr, char option){
    char jobvl = 'N';
    char jobvr = option;
    int n=m.rows();
    int lda=m.cols();
    assert(m.IsSquare());

    int ldvl = 1;
    int ldvr = n;
    std::vector<fcomplex> vl;
    evals.resize(n);
    vr.resize(ldvr * n);
    std::vector<fcomplex> work(3);
    std::vector<float> rwork(2 * n);
    int lwork= -1;
    int info;

    // query:
    LAPACK::cgeev_(&jobvl, &jobvr, &n, &(m(0,0)), &lda, &(evals[0]), &(vl[0]), &ldvl, &(vr[0]), &ldvr,
           &(work[0]), &lwork, &(rwork[0]), &info);
    if (info!=0) {
        std::cerr<<"info="<<info<<" during query\n";
        throw std::runtime_error("diag: cgeev_: failed with info!=0.\n");
    }
    const int NB = 256;
    lwork = std::max(1 + static_cast<int>(std::real(work[0])), (NB + 2)*n);
    work.resize(lwork);

    // real work:
    LAPACK::cgeev_(&jobvl, &jobvr, &n, &(m(0,0)), &lda, &(evals[0]), &(vl[0]), &ldvl, &(vr[0]), &ldvr,
           &(work[0]), &lwork, &(rwork[0]), &info);
    if (info!=0) {
        std::cerr<<"work info="<<info<<"\n";
        throw std::runtime_error("diag: cgeev_: failed with info!=0.\n");
    }
}


// ============================================================================
// =                           Matrix Scalar Wrapper                          =
// ============================================================================
void mscal(const dcomplex& a, Matrix<dcomplex>& X) {
    int incx = 1;
    int n = X.rows() * X.cols();
    BLAS::zscal_(&n, &a, &(X(0,0)), &incx);
}

void mscal(const fcomplex& a, Matrix<fcomplex>& X) {
    int incx = 1;
    int n = X.rows() * X.cols();
    BLAS::cscal_(&n, &a, &(X(0,0)), &incx);
}

void mscal(const double & a, Matrix<double>& X) {
    int incx = 1;
    int n = X.rows() * X.cols();
    BLAS::dscal_(&n, &a, &(X(0,0)), &incx);
}

void mscal(const float & a, Matrix<float>& X) {
    int incx = 1;
    int n = X.rows() * X.cols();
    BLAS::sscal_(&n, &a, &(X(0,0)), &incx);
}


// ============================================================================
// =                Matrix Vector Multiplication Wrapper                      =
// ============================================================================
// real matrix vector multiplication with single precision
void mxvw(Matrix<float>& A, std::vector<float>& X, std::vector<float>& Y, char option) {
    int incx = 1; int incy = 1; // stride
    if (option == 'g') {
        // general matrix
        char trans = 'N';
        float alpha = 1;
        float beta = 0;
        int m = A.rows();
        int n = A.cols();
        int lda = m;
        assert(X.size() == n);
        Y.resize(m);
        BLAS::sgemv_(&trans, &m, &n, &alpha, &(A(0,0)), &lda, &(X[0]), &incx, &beta, &(Y[0]), &incy);
    }
    else if (option == 's') {
        // symmetric matrix

    }
    else {
        std::cerr<<"option="<<option<<"\n";
        throw std::logic_error("diag: mxvw_float: invalid option.\n");
    }

}

// real matrix vector multiplication with double precision
void mxvw(Matrix<double>& A, std::vector<double>& X, std::vector<double>& Y, char option) {
    int incx = 1; int incy = 1; // stride
    if (option == 'g') {
        // general matrix
        char trans = 'N';
        double alpha = 1;
        double beta = 0;
        int m = A.rows();
        int n = A.cols();
        int lda = m;
        assert(X.size() == n);
        Y.resize(m);
        BLAS::dgemv_(&trans, &m, &n, &alpha, &(A(0,0)), &lda, &(X[0]), &incx, &beta, &(Y[0]), &incy);
    }
    else if (option == 's') {
        // symmetric matrix

    }
    else {
        std::cerr<<"option="<<option<<"\n";
        throw std::logic_error("diag: mxvw_double: invalid option.\n");
    }

}

// complex matrix vector multiplication with single precision
void mxvw(Matrix<fcomplex>& A, std::vector<fcomplex>& X, std::vector<fcomplex>& Y, char option) {
    int incx = 1; int incy = 1; // stride
    if (option == 'g') {
        // general matrix
        char trans = 'N';
        fcomplex alpha = 1;
        fcomplex beta = 0;
        int m = A.rows();
        int n = A.cols();
        int lda = m;
        assert(X.size() == n);
        Y.resize(m);
        BLAS::cgemv_(&trans, &m, &n, &alpha, &(A(0,0)), &lda, &(X[0]), &incx, &beta, &(Y[0]), &incy);
    }
    else if (option == 'h') {
        // hermitian matrix
        assert(A.IsHermitian());
        char uplo = 'U';
        int n = A.cols();
        int lda = n;
        Y.resize(n);
        fcomplex alpha = 1;
        fcomplex beta = 0;
        BLAS::chemv_(&uplo, &n, &alpha, &(A(0,0)), &lda, &(X[0]), &incx, &beta, &(Y[0]), &incy);
    }
    else {
        std::cerr<<"option="<<option<<"\n";
        throw std::logic_error("diag: mxvw_fcomplex: invalid option.\n");
    }

}

// complex matrix vector multiplication with double precision
void mxvw(Matrix<dcomplex>& A, std::vector<dcomplex>& X, std::vector<dcomplex>& Y, char option) {
    int incx = 1; int incy = 1; // stride
    if (option == 'g') {
        // general matrix
        char trans = 'N';
        dcomplex alpha = 1;
        dcomplex beta = 0;
        int m = A.rows();
        int n = A.cols();
        int lda = m;
        assert(X.size() == n);
        Y.resize(m);
        BLAS::zgemv_(&trans, &m, &n, &alpha, &(A(0,0)), &lda, &(X[0]), &incx, &beta, &(Y[0]), &incy);
    }
    else if (option == 'h') {
        // hermitian matrix
        assert(A.IsHermitian());
        char uplo = 'U';
        int n = A.cols();
        int lda = n;
        Y.resize(n);
        dcomplex alpha = 1;
        dcomplex beta = 0;
        BLAS::zhemv_(&uplo, &n, &alpha, &(A(0,0)), &lda, &(X[0]), &incx, &beta, &(Y[0]), &incy);

    }
    else {
        std::cerr<<"option="<<option<<"\n";
        throw std::logic_error("diag: mxvw_dcomplex: invalid option.\n");
    }

}


// ============================================================================
// =                Matrix Matrix Multiplication Wrapper                      =
// ============================================================================
// real matrix matrix multiplication with double precision
void mxmw(Matrix<double>& A, Matrix<double>& B, Matrix<double>& C, char option) {
    if (option == 'g') {
        char transa = 'N';
        char transb = 'N';
        int m = A.rows();  // rows of A, C
        int n = B.cols();  // cols of B, C
        C.resize(m, n);
        int k = A.cols(); // cols of A; rows of B
        double alpha = 1;
        double beta = 0;
        int lda = m;
        int ldb = k;
        int ldc = m;
        BLAS::dgemm_(&transa, &transb, &m, &n, &k, &alpha, &(A(0,0)), &lda, &(B(0,0)),
                     &ldb, &beta, &(C(0,0)), &ldc);
    }
    else {
        std::cerr<<"option="<<option<<"\n";
        throw std::logic_error("diag: mxmw_double: invalid option.\n");
    }
}

// real matrix matrix multiplication with single precision
void mxmw(Matrix<float>& A, Matrix<float>& B, Matrix<float>& C, char option) {
    if (option == 'g') {
        char transa = 'N';
        char transb = 'N';
        int m = A.rows();  // rows of A, C
        int n = B.cols();  // cols of B, C
        C.resize(m, n);
        int k = A.cols(); // cols of A; rows of B
        float alpha = 1;
        float beta = 0;
        int lda = m;
        int ldb = k;
        int ldc = m;
        BLAS::sgemm_(&transa, &transb, &m, &n, &k, &alpha, &(A(0,0)), &lda, &(B(0,0)),
                     &ldb, &beta, &(C(0,0)), &ldc);
    }
    else {
        std::cerr<<"option="<<option<<"\n";
        throw std::logic_error("diag: mxmw_float: invalid option.\n");
    }
}

// complex matrix matrix multiplication with double precision
void mxmw(Matrix<dcomplex>& A, Matrix<dcomplex>& B, Matrix<dcomplex>& C, char option) {
    if (option == 'g') {
        char transa = 'N';
        char transb = 'N';
        int m = A.rows();  // rows of A, C
        int n = B.cols();  // cols of B, C
        C.resize(m, n);
        int k = A.cols(); // cols of A; rows of B
        dcomplex alpha = 1;
        dcomplex beta = 0;
        int lda = m;
        int ldb = k;
        int ldc = m;
        BLAS::zgemm_(&transa, &transb, &m, &n, &k, &alpha, &(A(0,0)), &lda, &(B(0,0)),
                     &ldb, &beta, &(C(0,0)), &ldc);
    }
    else {
        std::cerr<<"option="<<option<<"\n";
        throw std::logic_error("diag: mxmw_dcomplex: invalid option.\n");
    }
}

// complex matrix matrix multiplication with single precision
void mxmw(Matrix<fcomplex>& A, Matrix<fcomplex>& B, Matrix<fcomplex>& C, char option) {
    if (option == 'g') {
        char transa = 'N';
        char transb = 'N';
        int m = A.rows();  // rows of A, C
        int n = B.cols();  // cols of B, C
        C.resize(m, n);
        int k = A.cols(); // cols of A; rows of B
        fcomplex alpha = 1;
        fcomplex beta = 0;
        int lda = m;
        int ldb = k;
        int ldc = m;
        BLAS::cgemm_(&transa, &transb, &m, &n, &k, &alpha, &(A(0,0)), &lda, &(B(0,0)),
                     &ldb, &beta, &(C(0,0)), &ldc);
    }
    else {
        std::cerr<<"option="<<option<<"\n";
        throw std::logic_error("diag: mxmw_fcomplex: invalid option.\n");
    }
}