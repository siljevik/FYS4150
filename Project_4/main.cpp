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
    int L = 20; // Lattize size
    int N = L*L; // Number of states/elements
    double E_empty = 0; // Initial energy
    double M_empty = 0; // Initial magnetism
    double J = 1.0; // Coupling constant = 1
    double beta =1/(T*k_b);
    // Creating an 'empty' matrix (filled with zeros)
    arma::mat S_empty(L, L);
    // Creating an empty list to fill with energies per atom from the lattice/matrix S
    vector<double> list_Es  = {};

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
    // We use this if we want an initial lattice with all spins up.
    arma::mat S = S_empty.fill(1); // Initially ORDERED
    // .. if not, we fill it up with random spins:
    //arma::mat S     = MCMC_s.spinnerboi(S_empty,L); // Initially UNORDERED

    // Calculating the total energy
    double E        = MCMC_s.tot_energyboi(S,L,E_empty,T,plusone,minusone,list_Es);
    // Calculating the total magnetism
    double M        = MCMC_s.tot_magnetboi(S,T,L,M_empty);


    /*============================================*/
    /*~~~~~      Analytical 2x2 lattice      ~~~~~*/
    /*============================================*/

    // Expected total energy J is the energy constant
    double Z        = analyticalboi.part_func(J, beta);

    // Expected total energy
    double exp_E    = analyticalboi.exp_tot_E(J,beta,Z); // Forventingsverdi for total energi
    double exp_EE   = analyticalboi.exp_tot_E_sqrd(J,beta,Z);
    // Expected total magentization
    double exp_M    = analyticalboi.exp_tot_M(J,beta,Z);
    double exp_MM   = analyticalboi.exp_tot_M_sqrd(J,beta,Z);
    // Specific heat capacity, CV, normalized to number of spins, N

    double CV       = analyticalboi.spec_heat_cap(N,k_b,T,exp_E,exp_EE);
    // Susceptibility, chi, normailzed to number of spins, N
    double chi      = analyticalboi.sus_chi(N,k_b,T,exp_M,exp_MM);

    //* ~~~~ Printing the analyticals ~~~~ *//
    /*cout << "\nAnalytical results: \n";
    cout << "\nExpected energy: " << exp_E;
    cout << "\n<epsilon>: "	<< exp_E/N;
    cout << "\n<|m|>: "		<< exp_M/N;
    cout << "\nCV: "		<< CV;
    cout << "\nX: "			<< chi << endl;*/
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//


    /*====================================*/
    /*~~~~~       Doing the MC       ~~~~~*/
    /*====================================*/

    // Initialize width and precision of the datafile
    int width = 15;                     /*====================================*//*====================================*//*====================================*/
    int prec = 6;                       /*====================================*//*====================================*//*====================================*/

    //open datafile to fill
    ofstream datafile;
    datafile.open("equilibrium_time_T_1_0.txt", ofstream::out | ofstream::trunc);
    //datafile.open("equilibrium_time_T_2_4.txt", ofstream::out | ofstream::trunc);
    //datafile.open("random_time_T_1_0.txt", ofstream::out | ofstream::trunc);
    //datafile.open("random_time_T_2_4.txt", ofstream::out | ofstream::trunc);

    //Header of txt-file
    datafile << "MC-cycles" << " " << "<epsilon>" << " " << "<|m|>" << endl; 

    double sum_E = 0; // Initial energy sum
    double sum_e = 0; // initial energy per spin sum
    int cycles = 8000;//pow(10,8); // Choosing how many MC cucles we want to do
    

    // Loop to count MCs
    for (int no = 0; no < cycles; no++)
    {
        // Cycle number
        int cycle = no + 1;
        // Empty vector for the new_E below
        vector<double> list_Es_new = {};

        arma::mat S_MC      = MCMC_s.random_spinnergal(S,T,L,N,E,M,beta,plusone,minusone);
        //cout << "\nEnergy before tot_energyboi" << E;
        
        double new_E 	    = MCMC_s.tot_energyboi(S_MC, L, E, T,plusone,minusone,list_Es_new);
        
        sum_E               += new_E;   //energy pr lattice (for avg)
        sum_e     	        += new_E/N; //energy pr spin (for avg)
        double E_avg        = sum_E/cycle; // mean energy///////////////////////////////
        double e_avg 	    = sum_e/cycle; // mean energy per spin
        cout << "\ne_avg: " << E_avg << " exp_E: " << exp_E;
        


        double new_M		= MCMC_s.tot_magnetboi_abs(S_MC, T, L, M);
        double exp_m	    = abs(new_M/N); // magnetization per spin
       
        


        //cout << "\nsum_E: " << sum_E << "\n";
        

        // Plot no som x og E_avg som y
        // now we print the results for the averge
        //cout << "mean energy(E): " << E_avg << " and cycle # " << no <<"\n";
        //cout << "mean energy(e): " << e_avg << "\n";
        //cout << "=======================" << endl;

        // Writing the values in the datafile
        datafile << cycle << " "<< e_avg << " " << exp_m << endl;

        }// end of for loop(no = cycles)


    datafile.close();




    // end of main
    return 0;
}