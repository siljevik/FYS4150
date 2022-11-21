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
#include <vector>
#include <map>

// including the header and function files
#include "MCMC_spin.hpp" // MCMC = Markov Chain Monte Carlo
#include "MCMC_spin.cpp"
#include "analytical.hpp"
#include "analytical.cpp"

using namespace std;


int main()
{
    //int argc, char* argv[]){    // Finnes eksempler på git i openMP parallellisering
    // Calling the analytical class
    MCMC_spin MCMC_s;
    analytical analyticalboi;

    /*====================================*/
    /*~~~~ Constants, Variables, etc. ~~~~*/
    /*====================================*/
    // Constants:
    double const k_b = 1;//Boltzman constant = 1, temperature has therefore energy dimension
    // Some known values
    // mysys.variabel = noe; altså fra class fil
    double T = 1.0; // Temperature   ----> Can be changed later
    int L = 2; // Lattize size
    int N = L*L; // Number of states/elements
    double E = 0; // Initial energy
    double M = 0; // Initial magnetism
    double J = 1.0; // Coupling constant = 1
    double beta =1/(T*k_b);
    // Creating an 'empty' matrix (filled with zeros)
    arma::mat S(L, L);

    // Creating a matrix with all up-spins
    arma::mat S_all_up(L,L);
    S_all_up.fill(1);
    ////////////
    // Creating an empty list to fill with energies per atom from the lattice/matrix S
    vector<double> list_Es  = {};
    /////////

    // Creating empty vectors to look up indexes so we can look
    // at neighbours of the states (in case of 'border'-states)
    vector<int> plusone{};
    vector<int> minusone{};
    // Then fill up plusone (e.g. [1,2,3, ... , L-1, 0])
    for (int i = 0; i<L-1; i++){plusone.push_back(i+1);}
    plusone.push_back(0);
    // Fill up minusone (e.g. [L-1,0,1,2, ... , L-2])
    minusone.push_back(L-1);
    for (int i = 0; i<L-1; i++){minusone.push_back(i);}



    /*=================================*/
    /*~~~~ Markov Chain Monte Calo ~~~~*/
    /*=================================*/
    // Filling the matrix up with random spins:
    arma::mat S2                            = MCMC_s.spinnerboi(S,L);
    // Calculating the total energy
    double E2                               = MCMC_s.tot_energyboi(S2,L,E,T,plusone,minusone,list_Es);
//    for(int k=0;k<N-1;k++){cout<<list_Es[k]<<"\n";} // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Creates a vector with energy for each atom
    //vector<double> tot_energy_pr_atom_list  =MCMC_s.energy_listboi(S2,L,plusone,minusone);
    // Calculating the total magnetism
    // double M2                               = MCMC_s.tot_magnetboi(S2,T,L,M);

    // The matrix we are doing calculations for (if it is very big we don't wanna print it)
/*
    if (L <= 10) {
        cout << "Matrix: \n" 		<< S2;}         // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    cout << "Total energy: " 	<< E2 << " J\n";    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
    // Printing Energylist (since it is a vector):
    //for(int i=0; i <tot_energy_pr_atom_list.size(); i++) { // To print a vector we must print one at the time
    //    cout <<tot_energy_pr_atom_list.at(i) <<' '; }



    /*=====================================*/
    /*~~~~~       Analytical 2x2      ~~~~~*/
    /*=====================================*/

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

    // Printing our expected energy and exp_M for a 2x2 lattice
    //cout << "Analytical results:"	<< "\n";
    //cout << "<epsilon>: "		<< exp_E/N 	<< "\n";
    //cout << "<|m|>: "			<< exp_M/N 	<< "\n";
    //cout << "CV: "			<< CV 		<< "\n";
    //cout << "X: "			<< chi 		<< endl;



    /*====================================*/
    /*~~~~~       Doing the MC       ~~~~~*/
    /*====================================*/

    // Initialize width and precision of the datafile
    int width = 15;
    int prec = 6;

    //open datafile to fill
    ofstream datafile;
    datafile.open("equilibrium_time_T_1_0.txt", ofstream::out | ofstream::trunc);

    datafile << setw(width) << "MC-cycles" << setw(width)
             << "<epsilon>" << setw(width)
             << "<|m|>"     << setw(width)
             << endl;

    double sum_E = 0; // Initial energy sum
    double sum_e = 0; // initial energy per spin sum

    int cycles = 15000; // Choosing how many MC cucles we want to do
    vector<int> cycle_nr = {};
    //double boltzman_n = 0; // Just making it 0
    // double boltzman_value = MCMC_s.boltzman_factors(beta,boltzman_n);
    //arma::mat S_MC(L,L);

    // Loop to count MCs
    for (int no = 0; no < cycles; no++)
    {
       	arma::mat S_MC      = MCMC_s.random_spinnergal(S_all_up,T,L,N,E,M,beta,plusone,minusone);
/*
        if (L <= 10)
        {
            cout << "Matrix: \n" << S_MC;           // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        }
*/
        vector<double> list_Es_new = {};

        double new_E 	    = MCMC_s.tot_energyboi(S_MC, L, E, T,plusone,minusone,list_Es_new);
        //cout << new_E;
   //     for(int i=0;i<N;i++){cout << list_Es_new[i] <<"\n";}// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        double new_exp_E    = new_E/N; // energy pr. spin
        //	cout << "new_E: " << new_E << "\n";
        sum_E               	+= new_E;
        sum_e     	        += new_E/N;
        double new_M		= MCMC_s.tot_magnetboi_abs(S_MC, T, L, M);
        double exp_m	    	= abs(new_M/N); // magnetization per spin
        //	cout <<"new_M: " << new_M << "\n";
        double E_avg        	= sum_E/cycles; // mean energy
        double e_avg 	    	= sum_e/cycles; // mean energy per spin
        //	cout << "Current sum E: " << sum_E << "\n";
        //    cout <<"Cycle number:" << no << " and avg. E:" << E_avg << endl;

        // Plot no som x og E_avg som y
        // now we print the results for the averge
        //cout << "mean energy(E): " << E_avg << " and " << no <<"\n";
        //cout << "mean energy(e): " << e_avg << "\n";
        //cout << "=======================" << endl;
        cycle_nr.push_back(no);
        // we want to write the values in our datafile
        datafile << setw(width) << setprecision(prec) << cycle_nr[no]
                << setw(width) << setprecision(prec) << e_avg
                << setw(width) << setprecision(prec) << exp_m
                << endl;
    }// end of for loop(no = cycles)
    datafile.close();

return 0;
}
