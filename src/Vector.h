/*
*
* Template Numerical Toolkit (TNT)
*
* Mathematical and Computational Sciences Division
* National Institute of Technology,
* Gaithersburg, MD USA
*
*
* This software was developed at the National Institute of Standards and
* Technology (NIST) by employees of the Federal Government in the course
* of their official duties. Pursuant to title 17 Section 105 of the
* United States Code, this software is not subject to copyright protection
* and is in the public domain. NIST assumes no responsibility whatsoever for
* its use by other parties, and makes no guarantees, expressed or implied,
* about its quality, reliability, or any other characteristic.
*
*/

//
// Adapted by shifeng on 3/13/21.
//

#ifndef CPPKIT_VECTOR_H
#define CPPKIT_VECTOR_H
#include<vector>
#include <complex>
#include <cassert>
#include <iostream>
#include <list>

typedef std::complex<double> dcomplex;
typedef std::complex<float> fcomplex;
typedef unsigned int Subscript;
template <class T>
class Vector {
public:
    typedef         T   value_type;
    typedef         T   element_type;
    typedef         T*  pointer;
    typedef         T*  iterator;
    typedef         T&  reference;
    typedef const   T*  const_iterator;
    typedef const   T&  const_reference;
    Subscript lbound() const { return 1;}

protected:
    T* v_;
    T* vm1_;        // pointer adjustment for optimzied 1-offset indexing
    Subscript n_;

    // internal helper function to create the array
    // of row pointers

    void initialize(Subscript N)
    {
        // adjust pointers so that they are 1-offset:
        // v_[] is the internal contiguous array, it is still 0-offset
        //
        assert(v_ == NULL);
        v_ = new T[N];
        assert(v_  != NULL);
        vm1_ = v_-1;
        n_ = N;
    }

    void copy(const T*  v)
    {
        Subscript N = n_;
        Subscript i;

#ifdef TNT_UNROLL_LOOPS
        Subscript Nmod4 = N & 3;
        Subscript N4 = N - Nmod4;

        for (i=0; i<N4; i+=4)
        {
            v_[i] = v[i];
            v_[i+1] = v[i+1];
            v_[i+2] = v[i+2];
            v_[i+3] = v[i+3];
        }

        for (i=N4; i< N; i++)
            v_[i] = v[i];
#else

        for (i=0; i< N; i++)
            v_[i] = v[i];
#endif
    }

    void set(const T& val)
    {
        Subscript N = n_;
        Subscript i;

#ifdef TNT_UNROLL_LOOPS
        Subscript Nmod4 = N & 3;
        Subscript N4 = N - Nmod4;

        for (i=0; i<N4; i+=4)
        {
            v_[i] = val;
            v_[i+1] = val;
            v_[i+2] = val;
            v_[i+3] = val;
        }

        for (i=N4; i< N; i++)
            v_[i] = val;
#else

        for (i=0; i< N; i++)
            v_[i] = val;

#endif
    }



    void destroy()
    {
        /* do nothing, if no memory has been previously allocated */
        if (v_ == NULL) return ;

        /* if we are here, then matrix was previously allocated */
        delete [] (v_);

        v_ = NULL;
        vm1_ = NULL;
    }


public:

    // access

    iterator begin() { return v_;}
    iterator end()   { return v_ + n_; }
    const iterator begin() const { return v_;}
    const iterator end() const  { return v_ + n_; }

    // destructor

    ~Vector()
    {
        destroy();
    }

    // constructors

    Vector() : v_(0), vm1_(0), n_(0)  {};

    Vector(const Vector<T> &A) : v_(0), vm1_(0), n_(0)
    {
        initialize(A.n_);
        copy(A.v_);
    }

    Vector(Subscript N, const T& value = T()) :  v_(0), vm1_(0), n_(0)
    {
        initialize(N);
        set(value);
    }

    Vector(Subscript N, const T* v) :  v_(0), vm1_(0), n_(0)
    {
        initialize(N);
        copy(v);
    }

    Vector(Subscript N, char *s) :  v_(0), vm1_(0), n_(0)
    {
        initialize(N);
        std::istringstream ins(s);

        Subscript i;

        for (i=0; i<N; i++)
            ins >> v_[i];
    }



    // methods
    //
    Vector<T>& resize(Subscript N)
    {
        if (n_ == N) return *this;

        destroy();
        initialize(N);

        return *this;
    }

    void print() {
        std::cout.precision(8);
        // std::cout << "shape:= (" << nrow << "," << ncol << ")" << std::endl;
        for(int i=0; i < n_; i++) {
            std::cout << v_[i] << "\t";
        }
    }

    void println() {
        std::cout.precision(8);
        // std::cout << "shape:= (" << nrow << "," << ncol << ")" << std::endl;
        for(int i=0; i < n_; i++) {
            std::cout << v_[i] << "\n";
        }
    }

    T max() {
        T value = v_[0];
        for (int i = 1; i < n_; ++i) {
            if (v_[i] > value) { value = v_[i]; }
        }
        return value;
    }

    void fill(T val) {
        for (int i = 0; i < n_; ++i) {
            v_[i] = val;
        }
    }

    void clear(){
        for(int i=0; i < n_; i++) {
            v_[i] = 0;
        }
    }

    // assignments
    //
    Vector<T>& operator=(const Vector<T> &A)
    {
        if (v_ == A.v_)
            return *this;

        if (n_ == A.n_)         // no need to re-alloc
            copy(A.v_);

        else
        {
            destroy();
            initialize(A.n_);
            copy(A.v_);
        }

        return *this;
    }

    Vector<T>& operator=(const T& scalar)
    {
        set(scalar);
        return *this;
    }

    inline Subscript dim() const
    {
        return  n_;
    }

    inline Subscript size() const
    {
        return  n_;
    }


    inline reference operator()(Subscript i)
    {
#ifdef TNT_BOUNDS_CHECK
        assert(1<=i);
        assert(i <= n_) ;
#endif
        return vm1_[i];
    }

    inline const_reference operator() (Subscript i) const
    {
#ifdef TNT_BOUNDS_CHECK
        assert(1<=i);
        assert(i <= n_) ;
#endif
        return vm1_[i];
    }

    inline reference operator[](Subscript i)
    {
#ifdef TNT_BOUNDS_CHECK
        assert(0<=i);
        assert(i < n_) ;
#endif
        return v_[i];
    }

    inline const_reference operator[](Subscript i) const
    {
#ifdef TNT_BOUNDS_CHECK
        assert(0<=i);
        assert(i < n_) ;
#endif
        return v_[i];
    }


};


/* ***************************  I/O  ********************************/

template <class T>
std::ostream& operator<<(std::ostream &s, const Vector<T> &A)
{
    unsigned int N=A.dim();

    s <<  N << "\n";

    for (unsigned int i=0; i<N; i++)
        s   << A[i] << " " << "\n";
    s << "\n";

    return s;
}

template <class T>
std::istream & operator>>(std::istream &s, Vector<T> &A)
{

    unsigned int N;

    s >> N;

    if ( !(N == A.size() ))
    {
        A.newsize(N);
    }


    for (unsigned int i=0; i<N; i++)
        s >>  A[i];


    return s;
}

// *******************[ basic matrix algorithms ]***************************


template <class T>
Vector<T> operator+(const Vector<T> &A,
                    const Vector<T> &B)
{
    Subscript N = A.dim();

    assert(N==B.dim());

    Vector<T> tmp(N);
    Subscript i;

    for (i=0; i<N; i++)
        tmp[i] = A[i] + B[i];

    return tmp;
}

template <class T>
Vector<T> operator-(const Vector<T> &A,
                    const Vector<T> &B)
{
    Subscript N = A.dim();

    assert(N==B.dim());

    Vector<T> tmp(N);
    Subscript i;

    for (i=0; i<N; i++)
        tmp[i] = A[i] - B[i];

    return tmp;
}

template <class T>
Vector<T> operator*(const Vector<T> &A,
                    const Vector<T> &B)
{
    Subscript N = A.dim();

    assert(N==B.dim());

    Vector<T> tmp(N);
    Subscript i;

    for (i=0; i<N; i++)
        tmp[i] = A[i] * B[i];

    return tmp;
}


template <class T>
T dot_prod(const Vector<T> &A, const Vector<T> &B)
{
    Subscript N = A.dim();
    assert(N == B.dim());

    Subscript i;
    T sum = 0;

    for (i=0; i<N; i++)
        sum += A[i] * B[i];

    return sum;
}



// ============================================================================
// =                          Declare Vector Scalar                           =
// ============================================================================
void vscal(const dcomplex& a, std::vector<dcomplex>& X);
void vscal(const fcomplex& a, std::vector<fcomplex>& X);
void vscal(const double& a, std::vector<double>& X);
void vscal(const float& a, std::vector<float>& X);

// ============================================================================
// =                           Declare Vector Norm                            =
// ============================================================================
double norm(const std::vector<dcomplex>& v);
float norm(const std::vector<fcomplex>& v);
double norm(const std::vector<double>& v);
float norm(const std::vector<float>& v);

// ============================================================================
// =                  Declare Vector Vector Multiplication                    =
// ============================================================================
// X dot Y
dcomplex dot(std::vector<dcomplex>& X, std::vector<dcomplex>& Y);
fcomplex dot(std::vector<fcomplex>& X, std::vector<fcomplex>& Y);
double dot(std::vector<double>& X, std::vector<double>& Y);
float dot(std::vector<float>& X, std::vector<float>& Y);

// X.conj dot Y
dcomplex cdot(std::vector<dcomplex>& X, std::vector<dcomplex>& Y);
fcomplex cdot(std::vector<fcomplex>& X, std::vector<fcomplex>& Y);


/*----------------End Declaration-----------------*/
// = Free functions: Vector rescale
template<class T>
std::vector<T> operator * (const std::vector<T>& X , const auto a) {
    T sa = static_cast<T>(a);  // from c++17
    std::vector<T> tmp(X);
    vscal(sa, tmp);
    return tmp;
}

template<class T>
void operator *= (std::vector<T>& X , const auto a) {
    T sa = static_cast<T>(a);
    vscal(sa, X);
}

// = Vector plus
template<class T>
std::vector<T> operator + (const std::vector<T>& X, const std::vector<T> &Y) {
    int n = X.size();
    assert(n == Y.size());
    std::vector<T> R(n);
    for (int i = 0; i < n; ++i) {
        R[i] = X[i] + Y[i];
    }
    return R;
}

template<class T>
void operator += (std::vector<T>& X, const std::vector<T> &Y) {
    int n = X.size();
    assert(n == Y.size());
    for (int i = 0; i < n; ++i) {
        X[i] = X[i] + Y[i];
    }
}

// = Vector minus
template<class T>
std::vector<T> operator - (const std::vector<T>& X, const std::vector<T> &Y) {
    int n = X.size();
    assert(n == Y.size());
    std::vector<T> R(n);
    for (int i = 0; i < n; ++i) {
        R[i] = X[i] - Y[i];
    }
    return R;
}

template<class T>
void operator -= (const std::vector<T>& X, const std::vector<T> &Y) {
    int n = X.size();
    assert(n == Y.size());
    for (int i = 0; i < n; ++i) {
        X[i] = X[i] - Y[i];
    }
}

#endif //CPPKIT_VECTOR_H
