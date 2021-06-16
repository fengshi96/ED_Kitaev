#ifndef CPPKIT_MATRIX_H
#define CPPKIT_MATRIX_H

#include <vector>
#include <cassert>
#include <iostream>
#include <random>
#include <typeinfo>
#include <iterator>
#include <functional>
#include <complex>
#include "BLAS.h"
#include "Vector.h"

typedef std::complex<double> dcomplex;
typedef std::complex<float> fcomplex;

template<class T>
class Matrix {
public:

    //set all elements to zero
    Matrix() : nrow(0), ncol(0)
    {}

    //allocate number of row col and elements
    Matrix(int nrow, int ncol) : nrow(nrow), ncol(ncol), data_(nrow * ncol)
    {}

    // initialize with matrix elements in col major
    Matrix(int nrow, int ncol, std::vector<T>& v) : nrow(nrow), ncol(ncol), data_(nrow * ncol)
    {data_ = v;}

    // copy constructor
    Matrix(const Matrix<T>& m) {
        nrow=m.nrow;
        ncol=m.ncol;
        data_=m.data_;
    }

    const T& operator()(int i,int j) const {
        //    std::cout << i << " " << j << std::endl;
#ifdef DEBUG
        assert(i < nrow && j < ncol);
        assert(i+ j * nrow < data_.size());
#endif
        return data_[i+ j * nrow];
    }

    T& operator()(int i, int j) {
        //    std::cout << i << " " << j << std::endl;
#ifdef DEBUG
        assert(i < nrow && j < ncol);
        assert(i+ j * nrow < data_.size());
#endif
        return data_[i+ j * nrow];
    }  // caller operator. returns a reference of that element

    void print() {
            std::cout.precision(8);
            //std::cout << "matrix shape:= (" << nrow << "," << ncol << ")" << std::endl;
            for(int i=0; i < nrow; i++) {
                for(int j=0; j < ncol; j++) {
                    std::cout << data_[i+ j * nrow] << "\t";
                }
                std::cout << std::endl;
            }
    }


    void clear(){
        for(int i=0; i < nrow; i++) {
            for(int j=0; j < ncol; j++) {
                data_[i+ j * nrow] = 0.0;
            }
        }
    }

    void resize(int newrow, int newcol) {
        assert(newrow>0 && newcol>0);
        nrow=newrow;
        ncol=newcol;
        data_.clear();
        data_.resize(newrow*newcol);
    }

    inline size_t size() {
        return data_.size();
    }
    inline int rows() { return nrow;}
    inline int cols() { return ncol;}

    void fillRand();
    void fill(T val) {
        std::fill(data_.begin(),data_.end(),val);
    }

    int numNonZeros(){
        int counter=0;
        for(int i=0; i<data_.size(); i++)
            if(data_[i]!=0.0) counter++;
        return counter;
    }

    void del(){
        nrow=0; ncol=0;
        data_.resize(0);
    }

    bool IsSquare();
    bool IsHermitian();
    bool IsSymmetric();

    void conj()
    {
        int n = data_.size();
        for (int i = 0; i < n; ++i)
            data_[i] = std::conj(data_[i]);
    }

    Vector<T> colspace(int col) {
        Vector<T> colvec(nrow);
        for (int r = 0; r < nrow; ++r) {
            colvec[r] = data_[r + col * nrow];
        }
        return colvec;
    }

    Matrix<T> operator + (const Matrix<T>& B);
    Matrix<T> operator - (const Matrix<T>& B);
    Matrix<T> operator * (T a);
    void operator += (const Matrix<T>& B);
    void operator -= (const Matrix<T>& B);
    void operator *= (T a);

    void transpose(Matrix<T>& m2, const Matrix<T>& m);
    void transpose();
    void ajoint(Matrix<T>& m2, const Matrix<T>& m);
    void ajoint();

    Matrix<T> prod(Matrix<T>& B, char option = 'g');
    std::vector<T> prod(std::vector<T>& X, char option = 'g');

    T trace();

private:
    int nrow, ncol;
    std::vector<T> data_;
};


// ############################################################################
// =                                 Functions                                =
// ############################################################################
template<class T>
Matrix<T> Matrix<T>::operator + (const Matrix<T>& B) {
    int m = B.nrow;
    int n = B.ncol;
    assert(this->nrow == m);
    assert(this->ncol == n);
    Matrix<T> tmp(m, n);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            tmp(i,j) = data_[i + j * nrow] + B(i,j);
        }
    }
    return tmp;
}

template<class T>
Matrix<T> Matrix<T>::operator - (const Matrix<T>& B) {
    int m = B.nrow;
    int n = B.ncol;
    assert(this->nrow == m);
    assert(this->ncol == n);
    Matrix<T> tmp(m, n);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            tmp(i,j) = data_[i + j * nrow] - B(i,j);
        }
    }
    return tmp;
}

template<class T>
Matrix<T> Matrix<T>::operator * (const T a) {
    Matrix<T> tmp(*this);
    mscal(a, tmp);
    return tmp;
}

template<class T>
void Matrix<T>::operator += (const Matrix<T>& B) {
    int m = B.nrow;
    int n = B.ncol;
    assert(this->nrow == m);
    assert(this->ncol == n);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            data_[i + j * nrow] = data_[i + j * nrow] + B(i,j);
        }
    }
}

template<class T>
void Matrix<T>::operator -= (const Matrix<T>& B) {
    int m = B.nrow;
    int n = B.ncol;
    assert(this->nrow == m);
    assert(this->ncol == n);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            data_[i + j * nrow] = data_[i + j * nrow] - B(i,j);
        }
    }
}

template<class T>
void Matrix<T>::operator *= (const T a) {
    mscal(a, *this);
}


template<class T>
bool Matrix<T>::IsSquare() {
    bool out = true;
    if (nrow != ncol) {
        out = false;
        return out;
    }
    return out;
}

template<class T>
bool Matrix<T>::IsHermitian() {
    double eps=1e-8;
    bool out=true;
    if (not IsSquare()) {out = false; return out;}
    for(int i=0; i < nrow; i++) {
        for(int j=0; j < i; j++) {
            T Hij = data_[i+ j * nrow];
            T Hji = data_[j+ i * nrow];
            if(std::norm(Hij - std::conj(Hji)) > eps) {
                std::string tmp = "Hij != Hji " + std::to_string(i)+"-"+ std::to_string(j)+" \n";
                std::cout << i << " \t " << j << " \t " << Hij << " \t " << Hji << std::endl;
                out=false;
                return out;
            }
        }
    }
    return out;
}

template<class T>
bool Matrix<T>::IsSymmetric() {
    double eps=1e-8;
    bool out = true;
    if (not IsSquare()) {out = false; return out;}
    for(int i=0; i < nrow; i++) {
        for(int j=0; j < i; j++) {
            T Hij = data_[i+ j * nrow];
            T Hji = data_[j+ i * nrow];
            if(Hij - Hji > eps) {
                std::string tmp = "Hij != Hji " + std::to_string(i)+"-"+ std::to_string(j)+" \n";
                std::cout << i << " \t " << j << " \t " << Hij << " \t " << Hji << std::endl;
                out = false;
                return out;
            }
        }
    }
    return out;
}

template<class T>
void Matrix<T>::transpose(Matrix<T>& m2,const Matrix<T>& m)
{
    m2.resize(m.cols(),m.rows());
    for (int i=0;i<m2.rows();++i)
        for (int j=0;j<m2.cols();++j)
            m2(i,j) = m(j,i);
}

template<class T>
void Matrix<T>::transpose() {
    if (IsSquare()) {
        for (int i = 0; i < nrow; ++i) {
            for (int j = 0; j < i; ++j) {
                T tmp = data_[i+ j * nrow];
                data_[i+ j * nrow] = data_[j+ i * nrow];
                data_[j+ i * nrow] = tmp;
            }
        }
    }
    else {
        int newnrow = ncol; int newncol = nrow;
        Matrix<T> m2(newnrow, newncol);
        for (int i=0; i<newnrow; ++i)
            for (int j=0; j<newncol; ++j)
                m2(i,j) = data_[j+ i * nrow];

        this->resize(newnrow, newncol);
        for (int i=0; i<newnrow; ++i)
            for (int j=0; j<newncol; ++j)
                data_[i+ j * nrow] = m2(i, j);
    }
}

template<class T>
void Matrix<T>::ajoint(Matrix<T>& m2, const Matrix<T>& m){
    m2(m);
    m2.conj();
    m2.transpose();
}

template<class T>
void Matrix<T>::ajoint(){
    this->conj();
    this->transpose();
}

template<class T>
T Matrix<T>::trace() {
    assert(IsSquare());
    T trace = 0.0;
    for (int i = 0; i < nrow; ++i){
        trace += data_[i+ i * nrow];
    }
    return trace;
}

template<class T>
void Matrix<T>::fillRand() {
    // First create an instance of an engine.
    std::random_device rnd_device;
    // Specify the engine and distribution.
    std::mt19937 mersenne_engine {rnd_device()};  // Generates random
    std::uniform_real_distribution<double> dist {0, 1};

    auto gen = [&dist, &mersenne_engine](){
        return dist(mersenne_engine);
    };
    generate(data_.begin(), data_.end(), gen);
} // https://stackoverflow.com/a/23143753/14853469


// = Matrix Vector Multiplication
template<class T>
std::vector<T> prod(Matrix<T>& A, std::vector<T>& X, char option = 'g') {
    // const std::type_info& ti = typeid(T);
    std::vector<T> Y;
    mxvw(A, X, Y, option);
    return Y;
}

// = Matrix Vector Multiplication v2
template<class T>
std::vector<T> Matrix<T>::prod(std::vector<T>& X, char option) {
    std::vector<T> Y;
    mxvw(*this, X, Y, option);
    return Y;
}

// = Matrix Matrix Multiplication
template<class T>
Matrix<T> prod(Matrix<T>& A, Matrix<T>& B, char option = 'g') {
    // const std::type_info& ti = typeid(T);
    Matrix<T> Y;
    mxmw(A, B, Y, option);
    return Y;
}

// = Matrix Matrix Multiplication v2
template<class T>
Matrix<T> Matrix<T>::prod(Matrix<T>& B, char option){
    Matrix<T> Y;
    mxmw(*this, B, Y, option);
    return Y;
}




// ############################################################################
// =   The following go to Matrix.cpp; The following go to Matrix.cpp;        =
// ############################################################################
// ============================================================================
// =                          Declare Matrix Scalar                           =
// ============================================================================
void mscal(const dcomplex& a, Matrix<dcomplex>& X);
void mscal(const fcomplex& a, Matrix<fcomplex>& X);
void mscal(const double& a, Matrix<double>& X);
void mscal(const float& a, Matrix<float>& X);


// ============================================================================
// =                  Declare Matrix Vector Multiplication                    =
// ============================================================================
void mxvw(Matrix<dcomplex>& A, std::vector<dcomplex>& X, std::vector<dcomplex>& Y, char option = 'g');
void mxvw(Matrix<fcomplex>& A, std::vector<fcomplex>& X, std::vector<fcomplex>& Y, char option = 'g');
void mxvw(Matrix<double>& A, std::vector<double>& X, std::vector<double>& Y, char option = 'g');
void mxvw(Matrix<float>& A, std::vector<float>& X, std::vector<float>& Y, char option = 'g');


// ============================================================================
// =                  Declare Matrix Matrix Multiplication                    =
// ============================================================================
void mxmw(Matrix<dcomplex>& A, Matrix<dcomplex>& B, Matrix<dcomplex>& C, char option = 'g');
void mxmw(Matrix<fcomplex>& A, Matrix<fcomplex>& B, Matrix<fcomplex>& C, char option = 'g');
void mxmw(Matrix<double>& A, Matrix<double>& B, Matrix<double>& C, char option = 'g');
void mxmw(Matrix<float>& A, Matrix<float>& B, Matrix<float>& C, char option = 'g');


// ============================================================================
// =                        Declare Diagonalization                           =
// ============================================================================
// diag real symmetric and Hermitian matrix
void diag(Matrix<dcomplex> &A, Vector<double> &evals, char option);
void diag(Matrix<fcomplex> &A, std::vector<float> &evals, char option);
void diag(Matrix<double> &A, std::vector<double> &evals, char option);
void diag(Matrix<float> &A, std::vector<float> &evals, char option);

// diag non-symmetric real and complex matrix
void diag(Matrix<dcomplex> &m, std::vector<dcomplex>& evals, std::vector<dcomplex>& vr, char option);
void diag(Matrix<fcomplex> &m, std::vector<fcomplex>& evals, std::vector<fcomplex>& vr, char option);
void diag(Matrix<double> &m, std::vector<double>& evalsRe, std::vector<double>& evalsIm,
          std::vector<double> vr, char option);
void diag(Matrix<float> &m, std::vector<float>& evalsRe, std::vector<float>& evalsIm,
          std::vector<float> vr, char option);

// ############################################################################
// =                        End of Declare for Matrix.cpp;                    =
// ############################################################################


#endif



