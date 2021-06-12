//
// Created by shifeng on 6/10/21.
//

/* ***************************************************
 *  Authors:Chris Bishop & Nirav D. Patel
 *	Data:January 14, 2016
 *	Platform: linux
 *  Language: C++
 *	Model: 3OrMC_MF-Test
 * ***************************************************
*/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <omp.h>
#include <vector>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <random>
#include <complex>
#include <cmath>
#include <cassert>
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>

using namespace std;
//#include "Vector.h"
#include "Memusage.h"
#include "Matrix.h"
#include "ParametersEngine.h"
#include "Lattice.h"
#include "QuantumBasis.h"
#include "Hamiltonian.h"
#include "Exactdiag.h"
#include "DynLanczosSolver.h"
#include "Observables.h"

typedef std::size_t SizeType;
typedef Eigen::VectorXf VectorXf;
typedef ConstVariables ConstVariablestype;
typedef QBasis Basistype;
typedef Hamiltonian HamilType;
typedef DynLanczosSolver DynLancType;
typedef Observables ObservType;

typedef std::complex<double> dcomp;    // your typedef
void PrintGSSz(Eigen::VectorXcd& Psi, Basistype& Basis, int& Nsite_);
void PrintGSConf(Eigen::VectorXcd& Psi, Basistype& Basis, int& Nsite_);
//void PrintGSConf(VectorXf& Psi, Basistype& Basis, int& Nsite_, int& basisLabel1);

int main(int argc, char *argv[]) {
    if (argc<2) {
        throw std::invalid_argument("USE:: executable inputfile");
    }

    string inputfile = argv[1];
    double pi=acos(-1.0);
    double vm, rss;
    dcomplex zero(0.0,0.0);
    dcomplex onei(0.0,1.0);

    //std::cout.unsetf ( std::ios::floatfield );
    //std::cout.precision(6);

    // Read from the inputfile
    ConstVariablestype Parameters(inputfile);
    omp_set_dynamic(0);
    omp_set_num_threads(Parameters.Threads);
    Lattice Lat(Parameters);
    int Nsite_ = Parameters.NumberofSites;


    // Build binary basis
    Basistype Basis(Parameters);
    int nn=Basis.basis.size();
    process_mem_usage(vm, rss);
    std::cout << "       VM (MB): " << int(vm/1024) << "       RSS (MB): " << int(rss/1024) << std::endl;


    // Build non-zero Diagonal and Tight-binding part of Hamiltonian
    HamilType Hamiltonian(Parameters,Lat,Basis);
    Hamiltonian.TightB_Ham();
    process_mem_usage(vm, rss);
    std::cout << "       VM (MB): " << int(vm/1024) << "       RSS (MB): " << int(rss/1024) << std::endl;

    Hamiltonian.Diagonal();
    process_mem_usage(vm, rss);
    std::cout << "       VM (MB): " << int(vm/1024) << "       RSS (MB): " << int(rss/1024) << std::endl;


    ObservType Observ(Parameters, Lat, Basis, Hamiltonian);
    DynLancType Lanczos(Parameters, Basis, Hamiltonian);
    Eigen::VectorXcd Psi;
    if(Parameters.Solver=="ED") {
        Psi = ExactDiag(Parameters, Basis, Hamiltonian, Observ);
    } else {
        Psi = Lanczos.Lanczos_Nirav(Parameters.LancType);
        cout << setprecision(6) << fixed;
        cout << " Ritz values: ";
        for(int i=0;i<Lanczos.TriDeval.size();i++) cout << Lanczos.TriDeval[i] << " ";
        cout << endl;
    }

    // -------------------------------
    // --- Ent Entropy ---
    // -------------------------------.
    cout << endl << endl;
    cout << "Calculating the topological entanglement entropy " << endl;


    // -- Define a cut;
    Eigen::VectorXi S1sites(6);  // typedef Matrix< int , Dynamic , 1> Eigen::VectorXi
    Eigen::VectorXi S2sites(6);
    S1sites << 0, 1, 2, 3, 4, 5;
    S2sites << 6, 7, 8, 9, 10, 11;
    //S1sites << 0, 7;
    //S2sites << 1, 2, 3, 4, 5, 6;
    Eigen::MatrixXi ABCutLabels = Basis.ABSplit(S1sites, S2sites);

    //cout << ABCutLabels << endl;
    cout << " --> Print ABCutLabels into ABCutLabels.dat";
    ofstream outfile;
    string outfilename1 = "ABCutLabels.dat";
    outfile.open(outfilename1);
    outfile << ABCutLabels;
    outfile.close();

    int S1Hil = pow(2,S1sites.size());
    int S2Hil = pow(2,S2sites.size());

    Eigen::MatrixXcd S1RDM(S1Hil,S1Hil);
    S1RDM.setZero();
    for (int i=0; i<S1Hil; i++) {
        for (int j=0; j<S1Hil; j++) {
            for (int k=0; k<S2Hil; k++) {
                SizeType ikL = ABCutLabels(i,k);
                SizeType jkL = ABCutLabels(j,k);
                S1RDM(i,j) += conj(Psi[ikL])*Psi[jkL];
            }
        }
    }
    cout << " --> Print reduced density matrix into S1RDM.dat\n";
    //cout << S1RDM << endl << endl;
    string outfilename2 = "S1RDM.dat";
    outfile.open(outfilename2);
    outfile << S1RDM;
    outfile.close();

    cout << " --> Trace of the reduced density matrix " << S1RDM.trace() << endl;



    Eigen::MatrixXcd tmp = S1RDM;
    Eigen::VectorXd D;
    Diagonalize('V',tmp,D);
    cout << " --> eigen values of the RDM = ";
    for(int i=0;i<D.size();i++) cout << D[i] << ", ";
    cout << endl;


    double entropy=0.0;
    for(int i=0;i<D.size();i++) {
        if(abs(D[i])>1e-6) {
            entropy += -D[i]*log(D[i]);
        }
    }
    cout << " --> vN Entropy = " << entropy << endl << endl;



    cout << "--------THE END--------" << endl;
}


/*=======================================================================
 * ======================================================================
 * ======================================================================
*/

void PrintGSSz(Eigen::VectorXcd& Psi, Basistype& Basis, int& Nsite_) {
    int k_maxup = Basis.basis.size();
    int n=k_maxup;

    Eigen::VectorXcd V(Nsite_+1);
    for(int i0=0;i0<n;i0++){
        int ket = Basis.basis[i0];
        dcomplex coef1 = Psi[ket];
        int SzT = Basis.NPartInState(ket,Nsite_);
        V[SzT] += coef1; //*conj(coef1);
    }
    cout << V << endl;
}


void PrintGSConf(Eigen::VectorXcd& Psi, Basistype& Basis, int& Nsite_) {
    int k_maxup = Basis.basis.size();
    int n=k_maxup;

    vector<pair<double,int> >V;
    for(int i=0;i<n;i++){
        pair<double,int>P=make_pair((Psi[i]*conj(Psi[i])).real(),i);
        V.push_back(P);
    }
    sort(V.begin(),V.end());

    for(int i0=n-1;i0>=n-10;i0--){
        int i = V[i0].second;
        int basisLabel1 = i; // i1*k_maxdn + i2;
        int ket = Basis.basis[i];
        double coef = V[i0].first; //Psi[basisLabel1];
        dcomplex coef1 = Psi[basisLabel1];

        cout << basisLabel1 << " - ";
        //if(coef*coef>0.02) {
        for(int s=0;s<Nsite_;s++){
            cout << " \\underline{";
            if(Basis.findBit(ket,s) == 1) {
                cout << " \\uparrow ";
            } else {
                cout << " \\downarrow ";
            }
            cout << "} \\ \\ ";
            if(s==Nsite_-1) cout << " \\ \\ ";
        }
        cout << "& \\ \\ ";
        cout << coef << " \\\\ " << coef1 << endl;
    }
    cout << "   " << endl;

}








