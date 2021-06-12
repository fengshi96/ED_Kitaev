//
// Created by shifeng on 5/11/21.
//
#include "Parser.h"

Parameters::Parameters(string inputfile)  {
    RandomSeed = matchstring(inputfile, "RandomSeed");
    RDMGeneratortype generator(RandomSeed);
    generator_ = generator;

    Model = matchstring2(inputfile, "Model");
    Geometry = matchstring2(inputfile, "Geometry");
    Solver = matchstring2(inputfile, "Solver");

    Lx = int(matchstring(inputfile, "NLx"));
    Ly = int(matchstring(inputfile, "NLy"));
    if (Model == "Kitaev") {
        NumberofSites = Lx * Ly * 2;
    }
    else {
        NumberofSites = Lx * Ly;
    }

//    Temperature = matchstring(inputfile, "Temperature");
//    LanczosEps = matchstring(inputfile, "LanczosEps");
    IsPeriodicX = bool(matchstring(inputfile, "IsPeriodicX"));
    IsPeriodicY = bool(matchstring(inputfile, "IsPeriodicY"));
//    LancType = matchstring2(inputfile, "LancType");
//    if (LancType == "OnlyEgs")
//        std::cout <<
//                  "WARNING:: Reortho will be turned off for LancType=OnlyEgs (set by default) " << std::endl;
    Threads = int(matchstring(inputfile, "Threads"));


//    // -- Lanczos Steps -- by default set to 200
//    try {
//        LanczosSteps = matchstring(inputfile, "LanczosSteps");
//    } catch (exception &e) {
//        LanczosSteps = 200;
//        std::cout << "LanczosSteps = " << LanczosSteps << " (set by default) " << std::endl;
//    }
//
//    try {
//        SpecialSite = matchstring(inputfile, "SpecialSite");
//    } catch (exception &e) {
//        SpecialSite = NumberofSites / 2;
//        std::cout << "SpecialSite = " << SpecialSite << " (set by default) " << std::endl;
//    }
}

double Parameters::matchstring(string file, const string match) {
    string test;
    string line;
    std::ifstream readFile(file);
    double amount;
    bool pass=false;

    if (!readFile.is_open())
        perror("error while opening file");

    while (std::getline(readFile, line)) {

        std::istringstream iss(line);
        if (std::getline(iss, test, '=')) {
            // ---------------------------------
            if (iss >> amount && test==match) {
                // cout << amount << endl;
                pass=true;
            }
            else {
                pass=false;
            }
            // ---------------------------------
            if(pass) break;
        }
    }

    if (readFile.bad())
        perror("error while reading file");

//    if (!pass) {
//        string errorout=match;
//        errorout+="= argument is missing in the input file!";
//        throw std::invalid_argument(errorout);
//    }
    std::cout << match << " = " << amount << std::endl;
    return amount;
} // ----------

string Parameters::matchstring2(string file, const string match) {
    string test;
    string line;
    std::ifstream readFile(file);
    string amount;
    bool pass=false;

    if (!readFile.is_open())
        perror("error while opening file");

    while (std::getline(readFile, line)) {
        std::istringstream iss(line);
        if (std::getline(iss, test, '=')) {
            // ---------------------------------
            if (iss >> amount && test==match) {
                // cout << amount << endl;
                pass=true;
            }
            else {
                pass=false;
            }
            // ---------------------------------
            if(pass) break;
        }
    }

    if (readFile.bad())
        perror("error while reading file");

//    if (!pass) {
//        string errorout=match;
//        errorout+="= argument is missing in the input file!";
//        throw std::invalid_argument(errorout);
//    }
    std::cout << match << " = " << amount << std::endl;
    return amount;
} // ----------

