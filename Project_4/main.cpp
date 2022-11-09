//////////////////////////////////
//
// Compile with: g++ -O2 main.cpp func/particle.cpp  -std=c++11 -I include -o main.exe -larmadillo
//                        ---> THIS ONE MUST BE EDITED <---
//
// Run with: ./main.exe
//
/////////////////////////////////
#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <chrono>
// including the header and function files containing our classes
#include "MCMC_spin.hpp"
#include "MCMC_spin.cpp"

using namespace std;


int main(){

    // Constants:
    double const k_b = 1.380649*pow(10.,-23.); // Boltzman constant [ m^2 * kg / (s^2 * K) ]

    // Known values
    double s_d = -1; // Spin up
    double s_u = 1; // Spin down
    double T = 1; // Temperature   ----> Can be changed later
    // Matrix with random spin up and downs
    double L = 2; // Lattize size
    arma::mat ising_grid = arma::mat(L, L);


return 0;
}