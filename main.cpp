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

int main(int argc, char *argv[]) {
    if (argc<2) {
        throw std::invalid_argument("USE:: executable inputfile");
    }
    std::string inputfile = argv[1];


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

    // Diagonalization
    Vector<dcomplex> Psi;
    Vector<double> Evals;
    if (Parameters.Solver == "ED") {
        diag(Hamil.HamMatrix, Evals, 'V');
        std::cout << "Beg of Eigen Value: " << std::endl;
        Evals.println();
        std::cout << "End of Eigen Value" << std::endl;
        Psi = Hamil.HamMatrix.colspace(0);  // ground state
    } else if (Parameters.Solver == "Lanczos") {

    } else {
      throw std::runtime_error ("Solver not recognized!");
    }


    // Entanglement
    // -- Define a cut;
    std::vector<SizeType> S1sites = {0, 1, 2};
    std::vector<SizeType> S2sites = {3, 4, 5, 6, 7};  // kept system

    Matrix<SizeType> ABCutLabels = Basis.ABSplit(S1sites, S2sites);

//    cout << endl << " --> Print ABCutLabels" << endl;
//    ABCutLabels.print();
//
//
    int S1Hil = pow(2,S1sites.size());
    int S2Hil = pow(2,S2sites.size());
//
    Matrix<dcomplex> S1RDM(S1Hil,S1Hil);
    for (int i=0; i<S1Hil; i++) {
        for (int j=0; j<S1Hil; j++) {
            for (int k=0; k<S2Hil; k++) {
                SizeType ikL = ABCutLabels(i,k);
                SizeType jkL = ABCutLabels(j,k);
                S1RDM(i,j) += conj(Psi[ikL])*Psi[jkL];
            }
        }
    }
    cout << "--> Printing RDM to S1RDM.dat\n";
    std::ofstream outfile;
    string outfilename2 = "S1RDM.dat";
    outfile.open(outfilename2);
    for (int i=0; i < S1RDM.size(); i++)
    {
        for (int j=0; j < S1RDM.size(); j++)
        {
            outfile << S1RDM(i,j);
        }
        outfile << "\n";
    }
    outfile.close();


     cout << "--> Trace of the reduced density matrix = " << S1RDM.trace() << endl;

     Vector<double> entSpec;
     diag(S1RDM, entSpec, 'V');
     cout << endl << "--> Beg of eigen values of RDM " << endl;
     entSpec.println();
     cout << "--> End of eigen values of RDM " << endl;

     double entropy=0.0;
     for(double i : entSpec) {
         if(i > 0) {
             entropy += -i*log(i);
         }
     }
     cout << "--> vN Entropy = " << entropy << endl << endl;



     cout << "--------THE END--------" << endl;


    return 0;
}
