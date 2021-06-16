//
// Created by shifeng on 6/14/21.
//

#ifndef BITHAMILTONIAN_LANCZOS_H
#define BITHAMILTONIAN_LANCZOS_H
Eigen::VectorXcd Lanczos_Nirav(string PsiOrAll){
    std::cout << "Begin Lanczos! --- "<<std::endl;
    int iter, Lexit;
    double E0, D0, diff;
    dcomplex zero(0.0,0.0);
    int STARTIT=3; //iteration which diagonz. begins
    int MAXiter=variables_.LanczosSteps;
    int EXITiter=-1;
    double LanczosEps=variables_.LanczosEps;
    double vm, rss;

    const int N=Basis_.basis.size();
    Eigen::VectorXcd V0(N), V1(N), V2(N), e(N), d(N);
    vector<dcomplex> alpha, beta;
    Matrix<dcomplex> TrEVec;
    vector<Eigen::VectorXcd> LancBasis(MAXiter+1);

    V0.setZero(); V1.setZero(); V2.setZero();
    e.setZero(); d.setZero();
    Vorig.resize(N); Vorig.setZero();
    Randomize(Vorig);			//Vorig = Vectpush_backorXf::Random(N);
    Vorig.normalize();			//normalize(Vorig);

//		if(N<512) {
//			return ExactDiag();
//		}

    // ----- Begin Lanczos -------------
    V0 = Vorig;
    //		for (int i=0; i<N; i++) {
    //			cout << V0[i] << endl;
    //		}
    //V0.fill(1.0);
    //V0.normalize();

    beta.push_back(zero);  //beta_0 not defined
    LancBasis[0] = V0;
    V1 = Hamil_.MatrixVectorMult(V0);		// |V1> = H |V0>
    dcomplex tmp = V0.dot(V1);
    alpha.push_back(tmp);

    V1 = V1 - alpha[0]*V0; // (Equation 2.3) Rev. Mod. Phys., 66, 3 (1994) - note beta = 0.0

    beta.push_back(V1.norm());
    V1.normalize(); //V1 = V1/beta[1];
    LancBasis[1] = V1;

    // done 0th iteration
    Lexit = 0;   //exit flag
    E0 = 1.0;    //previous iteration GS eigenvalue

    iter = 0;
//		// ------ while ----------------
    while(Lexit != 1) {

        iter++;
        V2 = Hamil_.MatrixVectorMult(V1); // V2 = H |V1>
        alpha.push_back(V1.dot(V2));
        V2 = V2 - alpha[iter]*V1 - beta[iter]*V0; // (Equation 2.3) Rev. Mod. Phys. , Vol. 66, No. 3, July 1994

        dcomplex norm = V2.norm();
        beta.push_back(norm);
        V2.normalize(); //V2 = V2/beta[iter+1];
        LancBasis[iter+1] = V2;

        // --- Reorthogonalization of basis states ---- //

        bool ReOrtho=true;
        if(PsiOrAll=="OnlyEgs") ReOrtho=false;
        if(ReOrtho && iter%1==0) {
            //cout << "  Re-orthogonalization " << endl;
            for(int i=0; i<iter+1; i++) {
                dcomplex rij = LancBasis[i].dot(LancBasis[iter+1]);

                // Vj = Vj - <Vi|Vj> Vi -- gram-schmid
#pragma omp parallel for
                for (int h=0;h<N;h++)
                    LancBasis[iter+1](h) -= LancBasis[i](h)*rij;
            }
            LancBasis[iter+1].normalize();
            V2 = LancBasis[iter+1];
        } //--------------------------------------------- //

        V0 = V1;
        V1 = V2;

        //diagonalize tri-di matrix
        D0 = TriDiag('N',alpha,beta,TrEVec, TriDeval); //TrEVec.print();
        diff = abs(E0-D0);

//            if(PsiOrAll=="OnlyEgs" && iter>2) {
//                LancBasis[iter-1].resize(1);
//            }
        process_mem_usage(vm, rss);

        std::cout << iter << " \t " <<setprecision(16)<< D0 <<  " \t "
                  << " Eps: " << diff
                  << " \t    VM (MB): " << int(vm/1024)
                  << "       RSS (MB): " << int(rss/1024)
                  << std::endl;

        if(diff < LanczosEps) {
            Lexit = 1;
            E0 = D0;
            std::cout << std::endl;
            std::cout << "  Output Energy: " <<setprecision(16)<< D0
                      << " After " << iter << " iterations "
                      << " Eps: " << diff
                      << std::endl;
            std::cout << std::endl;
        }
        else {
            E0=D0;
        }


        //			int pmax = (iter>6) ? 6: iter;
        //			for (int i=0; i<pmax; i++) {
        //				cout << TriDeval[i] << " ";
        //			}

        if(iter==MAXiter-1 || Lexit==1) {
            Lexit=1;
            D0 = TriDiag('V',alpha,beta,TrEVec, TriDeval); //TrEVec.print();
        }
    } // ------ while --------
    EXITiter = iter+1;

    std::cout << "  Lowest Eigen Value: " << D0 << endl;
    std::cout << "  With Lanczos EPS:   " << diff << endl;
    std::cout << "  Alpha and beta of Tridiag Matrix : ";
    for (int m=0; m<EXITiter; m++) {
        cout << "(" << real(alpha[m]) << "," << real(beta[m]) << ") ";
    }
    cout << endl << endl << flush;



    int loopval=100000000;
    Eigen::VectorXcd Psi(N);
    if(PsiOrAll=="OnlyEgs") {
        Psi.setZero();
        return Psi;
    } else if(PsiOrAll=="GS") {
        loopval=1;
        cout << "  creating all eigen-states " << endl;
        PsiAll.resize(loopval);
#pragma omp parallel for
        for (int k=0; k<loopval; k++) {
            PsiAll[k].resize(N);
            PsiAll[k].setZero();
            for (int h=0; h<N; h++) {
                for (int m=0; m<EXITiter; m++) {
                    dcomplex tmp = LancBasis[m](h)*TrEVec(m,k);
                    //dcomplex val(tmp.real(),0.0);
                    PsiAll[k](h) += tmp;
                    //cout << k << "-" << h << "-" << m << "-" << LancBasis[m](h) << "-" << TrEVec(m,k) << endl;
                }
            }
            process_mem_usage(vm, rss);
            cout << k << "    "
                 << " \t    VM (MB): " << int(vm/1024)
                 << "       RSS (MB): " << int(rss/1024)
                 << endl;
            cout.flush();
        }
        cout << endl;
        LancBasis.clear();
        Eigen::VectorXcd Psi = PsiAll[0];
        return Psi;
    } else if(PsiOrAll=="All") {
        loopval=EXITiter;
        cout << "  creating all eigen-states " << endl;
        PsiAll.resize(loopval);
#pragma omp parallel for
        for (int k=0; k<loopval; k++) {
            PsiAll[k].resize(N);
            PsiAll[k].setZero();
            for (int h=0; h<N; h++) {
                for (int m=0; m<EXITiter; m++) {
                    dcomplex tmp = LancBasis[m](h)*TrEVec(m,k);
                    //dcomplex val(tmp.real(),0.0);
                    PsiAll[k](h) += tmp;
                    //cout << k << "-" << h << "-" << m << "-" << LancBasis[m](h) << "-" << TrEVec(m,k) << endl;
                }
            }
            process_mem_usage(vm, rss);
            cout << k << "    "
                 << " \t    VM (MB): " << int(vm/1024)
                 << "       RSS (MB): " << int(rss/1024)
                 << endl;
            cout.flush();
        }
        cout << endl;
        LancBasis.clear();
        Eigen::VectorXcd Psi = PsiAll[0];
        return Psi;
    } else {
        throw runtime_error("DynLanczosSolver: Use PsiOrAll==GS or All");
    }

}
#endif //BITHAMILTONIAN_LANCZOS_H
