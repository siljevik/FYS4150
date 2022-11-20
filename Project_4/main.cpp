//////////////////////////////////////////////////////////////////////////////////////
//                                                                                  //
// Compile with: g++ -O2 main.cpp -std=c++11 -I include -o main.exe -larmadillo     //
//                                                                                  //
// Run with: ./main.exe                                                             //
//                                                                                  //
//////////////////////////////////////////////////////////////////////////////////////


// Bevause of int argc, char* argv[]
// ./main.exe  <L> <T_start> <T_end> <T_n_steps> <MC_cycles> <MC_cycles_burn_in>
//
// ex:
// ./main.exe  2 1.0 2.4

#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <random>
#include<vector>
#include <map>

// including the header and function files
#include "MCMC_spin.hpp" // MCMC = Markov Chain Monte Carlo
#include "MCMC_spin.cpp"
#include "analytical.hpp"
#include "analytical.cpp"

using namespace std;


int main(){ //int argc, char* argv[]){ // Finnes eksempler på git i openMP parallellisering

    // Calling the analytical class
    MCMC_spin MCMC_s;
//    analytical analyticalboi;

    /*=================================*/
    /*~~~~ Constants and Variables ~~~~*/
    /*=================================*/
    // Constants:
    double const k_b = 1;//Boltzman constant = 1, temperature has therefore energy dimension
    // Some known values
//    double s_d = -1; // Spin up
    // mysys.variabel = noe; altså fra class fil
//    double s_u = 1; // Spin down
    double T = 1; // Temperature   ----> Can be changed later
    int L = 2; // Lattize size
    double N = pow(L,2); // Number of states 
    double E = 0; // Initial energy
    double M = 0; // Initial magnetism
  //  double J = 1; // Coupling constant = 1
    double beta =1/(T*k_b);
    // Creating an 'empty' matrix (filled with zeros)
    arma::mat S(L, L);
    
    


    /*=================================*/
    /*~~~~ Markov Chain Monte Calo ~~~~*/
    /*=================================*/
    // Creating the plusone minusone vectors to be used
    //vector<int> plusone, minusone         = MCMC_s.plus_minus_bois(int L);
    // Filling the matrix up with random spins:
    arma::mat S2                            = MCMC_s.spinnerboi(S,L);
    // Calculating the total energy
    double E2                               = MCMC_s.tot_energyboi(S2,L,E,T);
    // Creates a vector with energy for each atom
    vector<double> tot_energy_pr_atom_list  = MCMC_s.energy_listboi(S2,L);
    // Calculating the total magnetism
//    double M2                               = MCMC_s.tot_magnetboi(S2,T,L,M);

    // The matrix we are doing calculations for (if it is very big we don't wanna print it)
    if (L <= 10) {
        cout << "Matrix: \n" << S2;}
    
    cout << "Total energy: " << E2 << " J\n";

    // Printing Energylist (since it is a vector):
    //for(int i=0; i <tot_energy_pr_atom_list.size(); i++) { // To print a vector we must print one at the time
    //    cout <<tot_energy_pr_atom_list.at(i) <<' '; }


    /*=====================================*/
    /*~~~~~       Analytical 2x2      ~~~~~*/
    /*=====================================*/
/*
    // Expected total energy J is the energy constant
    double Z        = analyticalboi.part_func(J, beta);

    // Expected total energy
    double exp_E    = analyticalboi.exp_tot_E(J,beta,Z);
    double exp_EE   = analyticalboi.exp_tot_E_sqrd(J,beta,Z);
    // Expected total magentization
    double exp_M    = analyticalboi.exp_tot_M(J,beta,Z);
    double exp_MM   = analyticalboi.exp_tot_M_sqrd(J,beta,Z);
    // Specific heat capacity, CV, normalized to number of spins, N
//    double CV       = analyticalboi.spec_heat_cap(N,k_b,T,exp_E,exp_EE);
    // Susceptibility, chi, normailzed to number of spins, N
//    double chi      = analyticalboi.sus_chi(N,k_b,T,exp_M,exp_MM);

    // Printing our expected energy for a 2x2 lattice
    cout << "exp_E: " << exp_E << "\n";
*/
    /*====================================*/
    /*~~~~~       Doing the MC       ~~~~~*/
    /*====================================*/
    double sum_E = 0; // For the energies
    int cycles = 10; // Choosing how many MC cucles we want to do
    //double boltzman_n = 0; // Just making it 0
    // double boltzman_value = MCMC_s.boltzman_factors(beta,boltzman_n);
    //arma::mat S_MC(L,L);

    // Loop to count MCs
    for (int no = 1; no < cycles; no++){
//        boltzman_value =
        arma::mat S_MC = MCMC_s.random_spinnergal(S,T,L,N,E,M,beta);
        sum_E += E;
        double E_avg = sum_E/no; // (i = MC_cycles)
        cout << no << " " << E_avg;
        // Plot i som x og E_avg som y
    }

    // Stuff me might need:
    /*
    arma::mat S_new = MCMC_s.random_spinnergal(S2,L);
    double E3 = MCMC_s.tot_energyboi(S_new,L,E2);
    double first_p_sT = MCMC_s.prob_func(beta,E2,E3,Z);
    double E_p_sT = first_p_sT*E2;
    double E_before = E3;
    
    


    //Checking if the p_sT is less than or else (equal to)
    while (E_p_sT > exp_E) // Må denne gjøres om pga vi aldri får 
    // nøyaktig lik exp_E
    {
        // We spin one random atom in the lattice
        arma::mat S_next = MCMC_s.random_spinnergal(S_new,L);
        // Renaming energy
        
        // Calculate the new energy
        double E_after = MCMC_s.tot_energyboi(S_next,L,E_before); // Setter inn E (=0) for å sjekke om det kan hjelpe
        // Calculate new p_sT
        double p_sT = MCMC_s.prob_func(beta,E_before,E_after,Z);
        // Calculating E_p_sT
        E_p_sT = p_sT*E_after; // E_before
        cout << "\n E_p_sT: " << E_p_sT;
        // Updating values
        MC_count += 1;
        S_new = S_next;
        E_before = E_after;

        
    }
    */

    
    

return 0;
}
