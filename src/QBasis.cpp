//
// Created by shifeng on 5/5/21.
//
#include <cassert>
#include "QBasis.h"
typedef std::size_t SizeType;


// Define class functions
void QBasis::Initialize() {
    if (parameters_.type_id == 0) {
        // Fermions
        nHil_ = pow(2,Nsite_);
        basis.resize(nHil_);
        for(int ket=0; ket<pow(2,Nsite_); ket++) {
            //PrintBitwise(ket);
            basis[ket] = ket;
        }
        std::cout << "number of ket basis: " << basis.size() << std::endl;
    }
    else if (parameters_.type_id == 1) {
        // Bosons
    }
    else if (parameters_.type_id == 2) {
        // SpinHalf
        nHil_ = pow(2,Nsite_);
        basis.resize(nHil_);
        for(int ket=0; ket<pow(2,Nsite_); ket++) {
            //PrintBitwise(ket);
            basis[ket] = ket;
        }
        std::cout << "number of ket basis: " << basis.size() << std::endl;
    }
    else {
        string errorout = "invalid Parameter_Type!";
        throw std::invalid_argument(errorout);
    }

}

void QBasis::ApplySz(const int& site, int& ketin, int& ketout, std::complex<double>& coef) {
    coef=std::complex<double>(0.0,0.0);
    ketout=-1;

    int bit = int(findBit(ketin,site));
    coef = (bit==0) ? std::complex<double>(0.5,0.0) : std::complex<double>(-0.5,0.0);;
    ketout=ketin;
    assert(ketout>-1);
}


void QBasis::ApplySx(const int& site, int& ketin, int& ketout, std::complex<double>& coef) {
    coef=std::complex<double>(0.5,0.0); //0.5+0i;
    ketout=-1;
    ketout = flipBit(ketin,site);
}

void QBasis::ApplySy(const int& site, int& ketin, int& ketout, std::complex<double>& coef) {
    coef=std::complex<double>(0,0); //0.+0.i;
    ketout=-1;
    int bit = int(findBit(ketin,site));
    ketout = flipBit(ketin,site);
    coef = (bit==0) ? std::complex<double>(0.0,0.5): std::complex<double>(0.0,-0.5); //0.+0.5i : 0.-0.5i;
}

Matrix<SizeType> QBasis::ABSplit(std::vector<SizeType>& S1sites, std::vector<SizeType>& S2sites) {
    SizeType S1 = pow(2,S1sites.size());
    SizeType S2 = pow(2,S2sites.size());
    assert(S1sites.size() + S2sites.size() == Nsite_);
    Matrix<SizeType> Labels(S1,S2);

    for(SizeType S1ket=0; S1ket<S1; S1ket++) {
        Vector<SizeType> binList(Nsite_);

        for(SizeType i=0; i<S1sites.size(); i++) {
            int site = S1sites[i];
            assert(site<Nsite_);
            binList[site] = findBit(S1ket, i);
        } // - close i

        for(SizeType S2ket=0; S2ket<S2; S2ket++) {

            for(SizeType i=0; i<S2sites.size(); i++) {
                int site = S2sites[i];
                assert(site<Nsite_);
                binList[site] = findBit(S2ket, i);
            } // - close i

            int fullHLabel = binVec2int(binList); // - must come after the i loop
            Labels(S1ket,S2ket) = fullHLabel;
        } // - close S2ket
    } // - close S1ket

    return Labels;
}

/*  ********************************************
    Here we create a mask, apply the mask to n,
    and then right shift the masked value to get just
    the bit we want. We could write it out more fully as:
    *
    int mask =  1 << position;  			create mask
    int masked_pos = i & position;			apply mask to i
    int thebit = masked_pos >> position;	right shift the masked value to get bit at position
*/
int findBit(const int i, const int position){
    return (i & ( 1 << position )) >> position;
} // ---------------

int flipBit(const int ketin, int position){
    return ketin^(1<<position);
} // ---------------

int hoppstate(int initial_state,int hopp_from,int hopp_to) {
    assert(findBit(initial_state,hopp_from)==1);
    assert(findBit(initial_state,hopp_to)==0);

    int mask = (1 << hopp_from) + (1 << hopp_to);
    return mask^initial_state;
}
