//////////////////////////////////
//
// Compile with: g++ -O2 main.cpp -std=c++11 -I include -o main.exe -larmadillo
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
#include <random>

// including the header and function files
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
    int L = 2; // Lattize size
    int N = pow(L,2); // number of states 
    // Creating an 'empty' matrix (filled with zeros)
    arma::mat S(L, L, fill::zeros);
    cout << S;


    // Calling the class
    MCMC_spin MCMC_s;
    // Filling the matrix up with random spins:
    MCMC_s.spinnerboi(S);

    MCMC_s.tot_energy(S,L); 

    cout << S;


/*
    MCMC_spin mysystem(2.4);


    mysystem.energy();

*/    
return 0;
}
