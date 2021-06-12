//
// Created by shifeng on 5/11/21.
//

#ifndef BITHAMILTONIAN_HAMILTONIAN_H
#define BITHAMILTONIAN_HAMILTONIAN_H

#include "QBasis.h"
#include "Lattice.h"
#include "Memusage.h"

class Hamiltonian {
public:
    Hamiltonian(SParameters &variables, Lattice &Lat, QBasis &Basis) :
            variables_(variables),
            Basis_(Basis),
            Lat_(Lat),
            Nsite_(variables.NumberofSites) {
        SetupConnectors();
        TightB_Ham(); process_mem_usage();
        Diagonal(); process_mem_usage();
        Build();
    }

    Matrix<dcomplex> HamMatrix;
    Matrix<int> Kxx_pair, Kyy_pair, Kzz_pair;
    std::vector<double> Kxx_Conn, Kyy_Conn, Kzz_Conn;
    std::vector<int> Sxx_init, Sxx_final;
    std::vector<int> Syy_init, Syy_final;
    std::vector<int> Szz_init, Szz_final;
    std::vector<dcomplex> HTSxx, HTSyy, HDiag;

    void SetupConnectors(){
        clear();
        if (variables_.Model=="Heisenberg") {
            Heisenberg_Connectors();
        } else if (variables_.Model=="Kitaev") {
            Kitaev_Connectors(); // Build Connectors in Single-particle states
        } else {
            std::cerr<<"Model="<<variables_.Model<<"\n";
            throw std::string("Unknown Model parameter \n");
        }
    }

    void clear() {
        Kxx_Conn.resize(0); Kyy_Conn.resize(0); Kzz_Conn.resize(0);
        Sxx_init.clear(); Syy_init.clear(); Szz_init.clear();
        Sxx_final.clear(); Syy_final.clear(); Szz_final.clear();
        HTSxx.clear(); HTSyy.clear(); HDiag.clear();
        Kxx_pair.clear(); Kyy_pair.clear(); Kzz_pair.clear();
    }

    void Kitaev_Connectors();
    void Heisenberg_Connectors();

    void TightB_Ham();
    void Diagonal();
    void Build();

private:
    SParameters& variables_;
    Lattice& Lat_;
    QBasis& Basis_;
    std::vector<double> Conn, TJH_Conn, Pair_Conn, TJHSz_Conn;
    int Nsite_,Nup_,Ndn_;
};





#endif //BITHAMILTONIAN_HAMILTONIAN_H
