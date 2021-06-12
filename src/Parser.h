//
// Created by shifeng on 5/5/21.
//

#ifndef BITHAMILTONIAN_PARSER_H
#define BITHAMILTONIAN_PARSER_H
#include <exception>
#include <iostream>
#include <string>
#include <random>
#include <fstream>
#include <sstream>

typedef std::string string;
typedef std::exception exception;
typedef std::mt19937_64 RDMGeneratortype;

struct Parameters {
    enum typeID {
        FERMION,
        BOSON,
        SPINHALF
    };

//    static Parameters *Create(int id);
    Parameters() = default;
    explicit Parameters(string inputfile);

    RDMGeneratortype generator_;
    std::uniform_real_distribution<double> dis_;

    int Lx{}, Ly{}, NumberofSites{}, SpecialSite{}, nHil{};
    int LanczosSteps{}, RandomSeed{}, Threads{};
    double Temperature{}, LanczosEps{};
    bool IsPeriodicX{}, IsPeriodicY{};
    string Geometry, Model, Solver, LancType;

//    bool LocalSkw,RIXSDynamics,SpinCurrDynamics,EnergyCurrDynamics;
//    bool FiniteTemp;
//    int RIXSIntermediate_Steps;

//    virtual void Initialize(string inputfile);
    double grnd() {
        return dis_(generator_);
    } // double random number generator

    static double matchstring(string file, const string match);
    static string matchstring2(string file, const string match);
};

// ======================== for spin models =================================
// ======================== for spin models =================================
// ======================== for spin models =================================
struct SParameters : Parameters {
    const int type_id = SPINHALF;
    SParameters() = default;
    explicit SParameters(string inputfile) : Parameters(inputfile) {
        std::cout << "____________________________________" << std::endl;
        std::cout << " - Reading the inputfile: " << inputfile << std::endl;
        std::cout << "____________________________________" << std::endl;

        Model = matchstring2(inputfile,"Model");
        Kxx = matchstring(inputfile,"Kxx");
        Kyy = matchstring(inputfile,"Kyy");
        Kzz = matchstring(inputfile,"Kzz");
        Bxx = matchstring(inputfile,"Bxx");
        Byy = matchstring(inputfile,"Byy");
        Bzz = matchstring(inputfile,"Bzz");
        nHil = pow(2, NumberofSites);
        std::cout << "____________________________________" << std::endl << std::endl;
    }
    double Kxx, Kyy, Kzz, Bxx, Byy, Bzz;

};

// ======================== for Fermion models =================================
// ======================== for Fermion models =================================
// ======================== for Fermion models =================================
struct FParameters : Parameters {
    const int type_id = FERMION;
    FParameters() = default;
    FParameters(const string& inputfile) {

    }

};



// ======================== for Fermion models =================================
// ======================== for Fermion models =================================
// ======================== for Fermion models =================================
struct BParameters : Parameters {
    const int type_id = BOSON;
    BParameters() = default;
    BParameters(const string& inputfile) {

    }

};

//BParameters::BParameters(const string& inputfile) : Parameters(inputfile) {
//
//}



/* ***********
*  Base Functions ------
*  ***********
*/

//Parameters *Parameters::Create(int id)
//{
//
//    if( id == 1 )
//    {
//        return new SParameters;
//    }
//    else if( id == 2 )
//    {
//        return new TBParameters;
//    }
//}





#endif //BITHAMILTONIAN_PARSER_H
