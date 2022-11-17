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
#include<vector>

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
    double const k_b = 1;//Boltzman constant = 1, temperature has therefore energy dimension

    // Known values
    double s_d = -1; // Spin up
    double s_u = 1; // Spin down
    double T = 1; // Temperature   ----> Can be changed later
    // Matrix with random spin up and downs
    int L = 2; // Lattize size
    double N = pow(L,2); // number of states 
    int E = 0;
    int M = 0;
    double J = 1; // Coupling constant = 1
    double beta =1/(T*k_b);
    
    // Creating an 'empty' matrix (filled with zeros)
    arma::mat S(L, L);
    


    /*=================================*/
    /*~~~~ Markov Chain Monte Calo ~~~~*/
    /*=================================*/

    // Filling the matrix up with random spins:
    arma::mat S2 = MCMC_s.spinnerboi(S,L);

    // Calculating the total energy
    double E2 = MCMC_s.tot_energyboi(S2,L,E);

    // Creates a vector with energy for each atom
    vector<double> tot_energy_pr_atom_list = MCMC_s.energy_listboi(S2,L);


    // Calculating the total magnetism
    double M2 = MCMC_s.tot_magnetboi(S2,L,M);

    // The matrix we are doing calculations for (if it is very big we don't wanna print it)
    if (L <= 10) {
        cout << "Matrix: \n" << S2;}
    
    cout << "Total energy: " << E2 << " J\n";
    /*
    cout << "Energylist:\n";
    for(int i=0; i <tot_energy_pr_atom_list.size(); i++) {
        cout <<tot_energy_pr_atom_list.at(i) <<' '; }
    cout << "\n";
    cout << "Magnetism: " << M2 << " unit\n";
    */

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
    double CV       = analyticalboi.spec_heat_cap(N,k_b,T,exp_E,exp_EE);
    // Susceptibility, chi, normailzed to number of spins, N
    double chi      = analyticalboi.sus_chi(N,k_b,T,exp_M,exp_MM);

    // Testing testing 1-2-3
    cout << "exp_E: " << exp_E << "\n";
    /*
    cout << "exp_M: " << exp_M << "\n";
    cout << "Critical temperature: " << CV << "\n";
    cout << "Chi: " << chi << "\n";
    */


    /////////////////////////////////////////
    // Running the single_spinnergal function once
    arma::mat S_new = MCMC_s.single_spinnergal(S2,L);
    double E3 = MCMC_s.tot_energyboi(S_new,L,E2);
    double first_p_sT = MCMC_s.prob_func(beta,E2,E3,Z);
    double E_p_sT = first_p_sT*E2;

    double E_before = E3;
    // Creating a variable to be used for counting Monte Carlo cycles
    int MC_count = 0;
    //Checking if the p_sT is less than or else (equal to)
    while (E_p_sT > exp_E) // Må denne gjøres om pga vi aldri får 
    // nøyaktig lik exp_E
    {
        // We spin one random atom in the lattice
        arma::mat S_next = MCMC_s.single_spinnergal(S_new,L);
        // Renaming energy
        
        // Calculate the new energy
        double E_after = MCMC_s.tot_energyboi(S_next,L,E_before);
        // Calculate new p_sT
        double p_sT = MCMC_s.prob_func(beta,E_before,E_after,Z);
        // Calculating E_p_sT
        E_p_sT = p_sT*E_before;
        // Updating values
        MC_count += 1;
        S_next = S_new;
        E_before = E_after;
    }
    // Men hva skjer etter?
    cout << "\n We did " << MC_count << " Monte Carlo cycles to get good agreement with the analytical result.\n";
    /////////////////////////////////////////
    

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