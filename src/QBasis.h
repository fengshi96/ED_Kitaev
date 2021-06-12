//
// Created by shifeng on 5/5/21.
//

#ifndef BITHAMILTONIAN_QBASIS_H
#define BITHAMILTONIAN_QBASIS_H
#include <complex>
#include <vector>
#include "Matrix.h"
#include "Vector.h"
#include "Parser.h"

typedef std::size_t SizeType;
typedef std::complex<double> dcomplex;
typedef std::complex<float> fcomplex;

class QBasis {
public:
    QBasis(SParameters& parameters) : parameters_(parameters), Nsite_(parameters.NumberofSites)
    {
        Initialize();
    }

    void Initialize();

    void ApplySx(const int& site, int& ketin, int& ketout, std::complex<double>& coef);
    void ApplySy(const int& site, int& ketin, int& ketout, std::complex<double>& coef);
    void ApplySz(const int& site, int& ketin, int& ketout, std::complex<double>& coef);

    inline int findBit(const int i, const int position){
        return (i & ( 1 << position )) >> position;
    }

    inline int flipBit(const int ketin, int position){
        return ketin^(1<<position);
    }

    inline int binVec2int(Vector<SizeType> inV) {
        SizeType length = inV.size();
        int out = 0;
        for(int i=0; i<length; i++) {
            int bin = inV[i];
            out += bin*pow(2,i);
        }
        return out;
    }

    void PrintBitwise(int i) {
        for(int j=0;j<Nsite_;j++){
            //std::cout << ((i & ( 1 << j )) >> j) << "  ";
            std::cout << findBit(i,j) << "  ";
        }
        std::cout << std::endl;
    }

    Matrix<SizeType> ABSplit(std::vector<SizeType>& S1sites, std::vector<SizeType>& S2sites);

    std::vector<int> basis;
    int nHil_;


private:
    SParameters& parameters_;
    int Nsite_;
};



#endif //BITHAMILTONIAN_QBASIS_H
