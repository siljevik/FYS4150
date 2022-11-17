//////////////////////////////////
//
// Compile with: g++ -O2 main.cpp -std=c++11 -I include -o main.exe -larmadillo
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
#include "analytical.hpp"
#include "analytical.cpp"

using namespace std;


int main(){
    
    // Calling the MCMC_spin class (MCMC = Markov Chain Monte Carlo)
    //MCMC_spin MCMC_s;
    // Calling the analytical class
    MCMC_spin MCMC_s;
    analytical analyticalboi;
    

    /*=================================*/
    /*~~~~ Constants and Variables ~~~~*/
    /*=================================*/
    // Constants:
    double const k_b = 1;//1.380649*pow(10.,-23.); // Boltzman constant [ m^2 * kg / (s^2 * K) ]

    // Known values
    double s_d = -1; // Spin up
    double s_u = 1; // Spin down
    double T = 1; // Temperature   ----> Can be changed later
    // Matrix with random spin up and downs
    int L = 2; // Lattize size
    int N = pow(L,2); // number of states 
    int E = 0;
    int M = 0;
    double J = 1;
    double beta =1/(T*k_b);
    //cout << "Beta is: " << beta << "\n";
    


    /*=================================*/
    /*~~~~ Markov Chain Monte Calo ~~~~*/
    /*=================================*/
    
    // Creating an 'empty' matrix (filled with zeros)
    arma::mat S(L, L);//, fill::zeros);

    // Filling the matrix up with random spins:
    arma::mat S2 = MCMC_s.spinnerboi(S,L);

    // Calculating the total energy
    double E2 = MCMC_s.tot_energyboi(S2,L,E);

    // Calculating the total magnetism
    double M2 = MCMC_s.tot_magnetboi(S2,L,M);

    // 
    cout << "Matrix: \n" << S2;
    cout << "Energy: " << E2 << " J\n";
    cout << "Magnetism: " << M2 << " unit\n";
    
    

    /*==================================*/
    /*~~~~~       Analytical       ~~~~~*/
    /*==================================*/
    
    // Expected total energy J is the energy constant
    double Z        = analyticalboi.part_func(J, beta);
    
    // Expected total energy
    double exp_E    = analyticalboi.exp_tot_E(J,beta,Z);
    double exp_EE   = analyticalboi.exp_tot_E_sqrd(J,beta,Z);
    // Expected total magentization
    double exp_M    = analyticalboi.exp_tot_M(J,beta,Z);
    double exp_MM   = analyticalboi.exp_tot_M_sqrd(J,beta,Z);
    // Specific heat capacity, CV, normalized to number of spins, N
    double CV       = analyticalboi.spec_heat_cap(N,J,beta,k_b,T,exp_E,exp_EE);
    // Susceptibility, chi, normailzed to number of spins, N
    double chi      = analyticalboi.sus_chi(N,J,beta,k_b,T,exp_M,exp_MM);
    
    // Testing testing 1-2-3
    cout << "Z: " << Z << "\n";
    cout << "exp_E: " << exp_E << "\n";
    cout << "exp_E: " << exp_EE << "\n";
    cout << "exp_M: " << exp_M << "\n";
    cout << "exp_MM: " << exp_MM << "\n";
    cout << "Critical temperature: " << CV << "\n";
    cout << "Chi: " << chi << "\n";
    

return 0;
}
/*
%%%%% NAVN PÅ IKKE-BINÆRE SØSKEN AV MAMMA/PAPPA 
(Istedet for onkel/tante):
- tankel
- onte
- pibling
- auncle
- pappas forvirret søsken
- mine forelsdres perverse søsken
- mammas problematiske søsken
- titti
*/