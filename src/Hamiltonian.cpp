//
// Created by shifeng on 5/11/21.
//

#include "Hamiltonian.h"

//----------------functions--------------------------
void Hamiltonian::Diagonal() {
    std::cout << "create Diagonal matrix elements" << std::endl;
    int Hsize = Basis_.basis.size();
    std::cout << "	hilbert-space size:  "<< Hsize << std::endl;

    // allocate diagonal term!
    HDiag.resize(Hsize);
    dcomplex izero=(0.0,0.0);
    std::fill(HDiag.begin(),HDiag.end(),izero);

    // Szi Sz_j |Psi> = distrop down and create up
    // (create down remove up)_i (create up remove down)_j

    for(int bond=0;bond<Kzz_Conn.size();bond++) {
        int si=Kzz_pair(bond,0);
        int sj=Kzz_pair(bond,1);
        if(Kzz_Conn.size()==1) break;
        assert(si!=sj);
#pragma omp parallel for
        for(int i1=0; i1<Hsize; i1++) {
            int ket = Basis_.basis[i1];
            int keto1,keto2;
            dcomplex Szi,Szj;

            Basis_.ApplySz(si,ket,keto1,Szi); assert(ket==keto1);
            Basis_.ApplySz(sj,keto1,keto2,Szj); assert(keto1==keto2);

            HDiag[i1] += Kzz_Conn[bond]*Szi*Szj;
        }
    }


    // --- ADD Magnetic Field in spin-z-direction
#pragma omp parallel for
    for(int i1=0; i1<Hsize; i1++) {
        int ket = Basis_.basis[i1];
        dcomplex SzT=0;

        for(int site=0;site<Nsite_;site++) {
            dcomplex Szi; int keto1;
            Basis_.ApplySz(site,ket,keto1,Szi); assert(ket==keto1);
            SzT+=Szi;
        }
        HDiag[i1] += variables_.Bzz*SzT;
    }
}


void Hamiltonian::TightB_Ham() {
    std::cout << "create Tight-Binding/Off-Diag Hamiltonian" << std::endl;
    int Hsize = Basis_.basis.size();
    std::cout << "	hilbert-space size:  "<< Hsize << std::endl;

    process_mem_usage();

    // allocate diagonal term!
    dcomplex izero=(0.0,0.0);
    std::fill(HTSxx.begin(),HTSxx.end(),izero);
    std::fill(HTSyy.begin(),HTSyy.end(),izero);

    process_mem_usage();

    // Sxi Sx_j |Psi> = distroy down and create up
    // (create down remove up)_i (create up remove down)_j
    for(int bond=0;bond<Kxx_Conn.size();bond++) {
        int si=Kxx_pair(bond,0);
        int sj=Kxx_pair(bond,1);
        if(Kxx_Conn.size()==1) break;
        assert(si!=sj);
        for(int i1=0; i1<Hsize; i1++) {
            int ket = Basis_.basis[i1];

            int keto1,keto2;
            std::complex<double> Sxi,Sxj;
            Basis_.ApplySx(si,ket,keto1,Sxi);
            Basis_.ApplySx(sj,keto1,keto2,Sxj);

            //cout << Sxi << " \t " << Sxj << " \t " << Sxi*Sxj << endl;
            Sxx_init.push_back(ket);
            Sxx_final.push_back(keto2);
            HTSxx.push_back(Sxi*Sxj*Kxx_Conn[bond]);

        }
    }

    process_mem_usage();


    // Syi Sy_j |Psi> = distrop down and create up
    // (create down remove up)_i (create up remove down)_j
    for(int bond=0;bond<Kyy_Conn.size();bond++) {
        int si=Kyy_pair(bond,0);
        int sj=Kyy_pair(bond,1);
        if(Kyy_Conn.size()==1) break;
        assert(si!=sj);
        for(int i1=0; i1<Hsize; i1++) {
            int ket = Basis_.basis[i1];

            int keto1,keto2;
            std::complex<double> Syi,Syj;
            Basis_.ApplySy(si,ket,keto1,Syi);
            Basis_.ApplySy(sj,keto1,keto2,Syj);

            //cout << Syi << " \t " << Syj << " \t " << Syi*Syj << endl;
            Syy_init.push_back(ket);
            Syy_final.push_back(keto2);
            HTSyy.push_back(Syi*Syj*Kyy_Conn[bond]);
        }
    }

    process_mem_usage();



    // --- ADD Magnetic Field in spin-x and spin-y direction
    // x-direction
    for(int i1=0; i1<Hsize; i1++) {
        int ket = Basis_.basis[i1];
        for(int site=0;site<Nsite_;site++) {
            int keto1=-1;
            //int bit = int(findBit(ketin,site));
            keto1 = Basis_.flipBit(ket,site);
            dcomplex coef=std::complex<double>(0.5,0);; //0.5+0i;
            Sxx_init.push_back(keto1);
            Sxx_final.push_back(ket);
            HTSxx.push_back(coef*variables_.Bxx);

        }
    }

    process_mem_usage();


    // y-direction
    for(int i1=0; i1<Hsize; i1++) {
        int ket = Basis_.basis[i1];
        for(int site=0;site<Nsite_;site++) {
            int keto1=-1;
            dcomplex coef=std::complex<double>(0,0); //0.+0.i;
            Basis_.ApplySy(site,ket,keto1,coef);
            dcomplex magy = variables_.Byy*coef;
            Syy_init.push_back(keto1);
            Syy_final.push_back(ket);
            HTSyy.push_back(magy);

        }
    }

    process_mem_usage();
}


void Hamiltonian::Kitaev_Connectors() {
    std::cout << "creating Connectors:" << std::endl;

    Matrix<double> SinglePart_Kxx, SinglePart_Kyy; // local
    Matrix<double> SinglePart_Kzz;
    SinglePart_Kxx.resize(Nsite_,Nsite_);
    SinglePart_Kyy.resize(Nsite_,Nsite_);
    SinglePart_Kzz.resize(Nsite_,Nsite_);

    for(int i=0;i<Nsite_;i++){ 	// ith site

        int j = Lat_.N1neigh_(i,0);
        if(i<j)  {
            int jx = Lat_.indx_(j);
            int jy = Lat_.indy_(j);
            assert(Lat_.Nc_(jx,jy)!=-1);
            SinglePart_Kzz(i,j) = variables_.Kzz;
        }

        j = Lat_.N1neigh_(i,2);
        if(i<j)  {
            int jx = Lat_.indx_(j);
            int jy = Lat_.indy_(j);
            assert(Lat_.Nc_(jx,jy)!=-1);
            SinglePart_Kxx(i,j) = variables_.Kxx;
        }

        j = Lat_.N1neigh_(i,1);
        if(i<j)  {
            int jx = Lat_.indx_(j);
            int jy = Lat_.indy_(j);
            assert(Lat_.Nc_(jx,jy)!=-1);
            SinglePart_Kyy(i,j) = variables_.Kyy;
        }
    }




    int Kxxbonds=SinglePart_Kxx.numNonZeros();
    int Kyybonds=SinglePart_Kyy.numNonZeros();
    int Kzzbonds=SinglePart_Kzz.numNonZeros();

    cout << " Single particle Kitaev - xx " << endl;
    SinglePart_Kxx.print();

    cout << " Single particle Kitaev - yy " << endl;
    SinglePart_Kyy.print();

    cout << " Single particle Kitaev - zz " << endl;
    SinglePart_Kzz.print();


    if(Kxxbonds==0) Kxxbonds=1;
    if(Kyybonds==0) Kyybonds=1;
    if(Kzzbonds==0) Kzzbonds=1;
    Kxx_pair.resize(Kxxbonds,2); Kyy_pair.resize(Kyybonds,2); Kzz_pair.resize(Kzzbonds,2);
    Kxx_Conn.resize(Kxxbonds);   Kyy_Conn.resize(Kyybonds);   Kzz_Conn.resize(Kzzbonds);

    int counter=0;
    for(int i=0;i<Nsite_;i++) {
        for(int j=0;j<Nsite_;j++) {
            if(SinglePart_Kxx(i,j)!=0.0) {
                //assert(j>i);
                Kxx_pair(counter,0) = i;
                Kxx_pair(counter,1) = j;
                Kxx_Conn[counter] = SinglePart_Kxx(i,j);
                counter++;
            }
        }
    }

    counter=0;
    for(int i=0;i<Nsite_;i++) {
        for(int j=0;j<Nsite_;j++) {
            if(SinglePart_Kyy(i,j)!=0.0) {
                //assert(j>i);
                Kyy_pair(counter,0) = i;
                Kyy_pair(counter,1) = j;
                Kyy_Conn[counter] = SinglePart_Kyy(i,j);
                counter++;
            }
        }
    }

    counter=0;
    for(int i=0;i<Nsite_;i++) {
        for(int j=0;j<Nsite_;j++) {
            if(SinglePart_Kzz(i,j)!=0.0) {
                //assert(j>i);
                Kzz_pair(counter,0) = i;
                Kzz_pair(counter,1) = j;
                Kzz_Conn[counter] = SinglePart_Kzz(i,j);
                counter++;
            }
        }
    }
}

void Hamiltonian::Build() {
    int N = Basis_.basis.size();
    HamMatrix.resize(N, N);

    //make the Hamiltonian
    for (int i=0; i<N; i++) {
        dcomplex Hij=HDiag[i];
        HamMatrix(i,i)+=Hij;
    }

    int hilbert_t=HTSxx.size();
    // Kxx Kitaev - Sector -------
    for(int i=0;i<hilbert_t;i++){
        int Hi=Sxx_init[i];
        int Hj=Sxx_final[i];
        assert(Hi<N && Hj<N);
        dcomplex Hij=HTSxx[i];
        HamMatrix(Hi,Hj) += Hij;
    }

    hilbert_t=HTSyy.size();
    // Kxx Kitaev - Sector -------
    for(int i=0;i<hilbert_t;i++){
        int Hi=Syy_init[i];
        int Hj=Syy_final[i];
        assert(Hi<N && Hj<N);
        dcomplex Hij=HTSyy[i];
        HamMatrix(Hi,Hj) += Hij;
    }
}


void Hamiltonian::Heisenberg_Connectors() {
    std::cout << "FIXME11" << std::endl;
}