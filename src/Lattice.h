//
// Created by shifeng on 5/11/21.
// All credit go to Nirav D. Patel. All error goes to Shi Feng
//

#ifndef BITHAMILTONIAN_LATTICE_H
#define BITHAMILTONIAN_LATTICE_H

#include <cassert>
#include <iostream>
#include "Matrix.h"
#include "Vector.h"
#include "Parser.h"

using std::cout;
using std::endl;


class Lattice {
public:
    Lattice(Parameters& variables)
            : variables_(variables)
    {
        Initialize();
    }

    Matrix<int> Nc_,N1neigh_,N2neigh_;
    Vector<int> indx_,indy_;

    /*
     * ***********
     *  Functions in Class Coordinates ------
     *  ***********
    */
    void Initialize(){

        if(variables_.Geometry=="Honeycomb"){
            cout << "[Lattice.h]: Initializing Honeycomb lattice..." << endl;
            Nsite_ = variables_.NumberofSites;
            Number1Neigh_=3;
            Number2Neigh_=3;
            N1neigh_.resize(Nsite_,Number1Neigh_);
            N2neigh_.resize(Nsite_,Number2Neigh_);
            N1neigh_.fill(-1);
            N2neigh_.fill(-2);
            Honeycomb1();

        } else if(variables_.Geometry=="Chain"){
            Nsite_ = variables_.NumberofSites;
            Number1Neigh_=2;
            Number2Neigh_=2;
            N1neigh_.resize(Nsite_,Number1Neigh_);
            N2neigh_.resize(Nsite_,Number2Neigh_);
            N1neigh_.fill(-1);
            N2neigh_.fill(-2);
            //Chain();
        }
        else if (variables_.Geometry=="Square") {
            Nsite_ = variables_.NumberofSites;
            Number1Neigh_=4;
            if (variables_.IsPeriodicX || variables_.IsPeriodicY)
                Number1Neigh_=2;
            N1neigh_.resize(Nsite_,Number1Neigh_);
            N1neigh_.fill(-1);
            //Square();
        }
    }



    // ---------------- chain ---------------------------



    // ---------------- Honeycomb ---------------------------
    void Honeycomb1() {


        int LLX = variables_.Lx;
        int LLY = variables_.Ly;
        string shift = "right"; // 'right'
        std::cout << "creating Honeycomb lattice: " << shift << " shifted " << std::endl;
        double scalex = 2.0, scaley = 4.0 / sqrt(3.0);
        Vector<double> t1(2), t2(2);
        t1[0] = 1.0 * scalex;
        t1[1] = 0;

        if (shift == "left") {
            t2[0] = -0.5 * scalex;
            t2[1] = sqrt(3.0) / 2.0 * scaley;
        } else if (shift == "right") {
            t2[0] = 0.5 * scalex;
            t2[1] = sqrt(3.0) / 2.0 * scaley;
        }

        // Site labeling
        indx_.resize(Nsite_);
        indy_.resize(Nsite_);
        Nc_.resize(LLX * 2 + LLY, LLY * 2);
        Nc_.fill(-1);

        int xv = 0, counter = 0;
        for (int i = 0; i < LLX; i++) {
            if (i != 0) { xv += t1[0]; }
            int y0 = 0;
            int x1, y1;
            if (shift == "right") {
                x1 = xv + 1.0;
                y1 = 1.0;
            } else if (shift == "left") {
                x1 = xv - 1.0;
                y1 = 1.0;
            }

            for (int j = 0; j < LLY; j++) {
                int cxa = int(xv + j * t2[0]);
                int cxb = int(x1 + j * t2[0]);
                int cya = int(y0 + j * t2[1]);
                int cyb = int(y1 + j * t2[1]);

                indx_[counter] = cxa;
                indy_[counter] = cya;
                Nc_(cxa, cya) = counter;
                counter++;

                indx_[counter] = cxb;
                indy_[counter] = cyb;
                Nc_(cxb, cyb) = counter;
                counter++;
            }
        }
        Nc_.print();


        int xmax = indx_.max();
        int ymax = indy_.max();
        int jx, jy;
        for (int i = 0; i < Nsite_; i++) {    // ith site
            int ix = indx_[i];
            int iy = indy_[i];

            // SxSx (kitaev) - neighbor 0
            jx = ix + 1;
            jy = iy + 1;
            if (jx <= xmax && jy <= ymax && Nc_(jx, jy) != -1) {
                int j = Nc_(jx, jy);
                N1neigh_(i, 2) = j;
                N1neigh_(j, 2) = i;
            }

            // SySy (kitaev) - neighbor 1
            jx = ix + 1;
            jy = iy - 1;
            if(jx<=xmax && jy<=ymax && jy>=0 && Nc_(jx,jy)!=-1)  {
                int j = Nc_(jx,jy);
                N1neigh_(i,1) = j;
                N1neigh_(j,1) = i;
            }

            // SzSz (kitaev) - neighbor 2
            jx = ix;
            jy = iy + 1;
            if(jx<=xmax && jy<=ymax && Nc_(jx,jy)!=-1)  {
                int j = Nc_(jx,jy);
                N1neigh_(i,0) = j;
                N1neigh_(j,0) = i;
            }



            if(variables_.IsPeriodicY==true && shift=="right") {
                jx = int(ix - LLY);
                jy = 0;
                if(jx>=0 && iy==ymax && Nc_(jx,jy)!=-1) {
                    int j = Nc_(jx,jy);
                    N1neigh_(i,0) = j;
                    N1neigh_(j,0) = i;
                }
            } else if(variables_.IsPeriodicY==true && shift=="left") {
                jx = int(ix + LLY);
                jy = 0;
                if(jx<xmax+1 && iy==ymax and Nc_(jx,jy)!=-1) {
                    int j = Nc_(jx,jy);
                    N1neigh_(i,0) = j;
                    N1neigh_(j,0) = i;
                }
            }

            if(variables_.IsPeriodicX==true && shift=="right" && iy%2==0 && Nc_(ix,iy)<LLY*2) {
                jx = int(ix + LLX*2 - 1);
                jy = int(iy + 1);
                if(jx<=xmax && iy<=ymax && iy % 2 == 0 && Nc_(jx,jy)!=-1) {
                    int j = Nc_(jx,jy);
                    N1neigh_(i,1) = j;
                    N1neigh_(j,1) = i;
                }
            } else if(variables_.IsPeriodicX==true && shift=="left" && iy%2==1 && Nc_(ix,iy)<LLY*2) {
                jx = int(ix + LLX*2 - 1);
                jy = int(iy - 1);
                if(jx<xmax+1 && iy<ymax+1 && Nc_(jx,jy)!=-1) {
                    int j = Nc_(jx,jy);
                    N1neigh_(i,2) = j;
                    N1neigh_(j,2) = i;
                }
            }
        }


        cout << " 1st Nearest neighbors " << endl;
        N1neigh_.print();
        cout << endl;

//		cout << " 2nd Nearest neighbors " << endl;
//		N2neigh_.print();
//		cout << endl;

    } // end Honeycomb1



private:
    Parameters& variables_;
    int Nsite_,Number1Neigh_,Number2Neigh_;
};





#endif //BITHAMILTONIAN_LATTICE_H
