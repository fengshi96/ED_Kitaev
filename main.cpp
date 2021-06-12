#include <iostream>
#include <fstream>
#include <omp.h>
#include "src/Parser.h"
#include "src/Lattice.h"
#include "src/Vector.h"
#include "src/QBasis.h"
#include "src/Memusage.h"
#include "src/Hamiltonian.h"
using std::cout;
using std::fstream;
using std::endl;
typedef std::size_t SizeType;
double pi=acos(-1.0);

int main() {
    string inputfile = "/Barn/Lab/BitHamiltonian/input.inp";


    //std::cout.unsetf ( std::ios::floatfield );
    //std::cout.precision(6);

    // Read from the inputfile
    SParameters Parameters(inputfile);
    omp_set_dynamic(0);
    omp_set_num_threads(Parameters.Threads);
    Lattice Lat(Parameters);


    // Build binary basis
    QBasis Basis(Parameters);
    process_mem_usage();

    // Build non-zero Diagonal and Off Diagonal part of Hamiltonian
    Hamiltonian Hamil(Parameters, Lat, Basis);

    // Exact diagonalization
    Vector<double> Evals;
    diag(Hamil.HamMatrix, Evals, 'V'); Evals.println();
    Vector<dcomplex> Psi(Hamil.HamMatrix.colspace(0));  // ground state



    // Entanglement !!!! This doesn't work! I'm going out and have fun!
    // -- Define a cut;
//    std::vector<SizeType> S1sites = {0, 1, 2};
//    std::vector<SizeType>S2sites = {3, 4, 5, 6, 7};
//    Matrix<SizeType> ABCutLabels;
//    Basis.ABSplit(S1sites, S2sites);
//;
//    //cout << ABCutLabels << endl;
//    cout << endl << " --> Print ABCutLabels" << endl;
//    ABCutLabels.print();
//
//
//    int S1Hil = pow(2,S1sites.size());
//    int S2Hil = pow(2,S2sites.size());
//
//    Matrix<dcomplex> S1RDM(S1Hil,S1Hil);
//    for (int i=0; i<S1Hil; i++) {
//        for (int j=0; j<S1Hil; j++) {
//            for (int k=0; k<S2Hil; k++) {
//                SizeType ikL = ABCutLabels(i,k);
//                SizeType jkL = ABCutLabels(j,k);
//                S1RDM(i,j) += conj(Psi[ikL])*Psi[jkL];
//            }
//        }
//    }
//    cout << " --> Print reduced density matrix into S1RDM.dat\n";
//    S1RDM.print();


    //cout << S1RDM << endl << endl;
//    string outfilename2 = "S1RDM.dat";
//    outfile.open(outfilename2);
//    outfile << S1RDM;
//    outfile.close();

    // cout << " --> Trace of the reduced density matrix " << S1RDM.trace() << endl;




//    SParameters Sp(filename);
//    std::cout << Sp.type_id << std::endl;

//    Lattice lattice(Sp);

//    fstream readFile(filename);
//    if (!readFile.is_open())
//        perror("error while opening file");
//
//    string line;
//    while (getline(readFile, line)) {
//        cout << "Success!" << endl;
//    }

    return 0;
}
