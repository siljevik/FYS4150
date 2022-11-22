///////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                           //
// Compile with: g++ -O2 main.cpp -fopenmp -std=c++11 -I include -o main.exe -larmadillo     //
//                                                                                           //
// Run with: ./main.exe <T_min> <T_max> <n_T> <n_cycles> <Lattice> <output filename>         //
//                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////

#include "omp.h"  // OpenMP header
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <assert.h> //???
#include <cmath>
#include <map>
#include <vector>
#include <chrono>

using namespace std;

// Here we want to parallelize with the outer loop over temperature values
// (<T_min> <T_max>) over n_T times for temp. variable T


// the source of inspiration for the outer loop can be found hier:
// https://github.com/anderkve/FYS3150/blob/master/code_examples/omp_parallelization/main_omp_outer_loop.cpp

int main(int argc, const char* argv[]){

    // check number of command line arguments
    assert(argc == 7);

    // Read command line arguments
    const double T_min      = atof(argv[1]); 
    const double T_max      = atof(argv[2]);
    const int    n_T        = atoi(argv[3];) 
    const int    n_cycles   = atoi(argv[4]);
    const int    Lattice    = atoi(argv[5]);
    const string output_file_name = argv[5];

    //
    // Outer loops over values of T
    //

    const double delta_T = (T_max - T_min)/(n_T - 1); // n_T points corresponds to  (n_T -1) intervals


    #pragma omp parallel // start of parallel region
    {
        // prepare file
        ofstream datafile;
        // Giving each thread a filename
        const int my_thread = omp_get_thread_num();
        string my_file = output_file_name + ".thread_" + to_string(my_thread);
        datafile.open(my_file.c_str(), ofstream::out | ofstream::trunc);
        
        /*====================================*/
        /*~~~~ Constants, Variables, etc. ~~~~*/
        /*====================================*/
        // Constants:
        double const k_b = 1;//Boltzman constant = 1, temperature has therefore energy dimension
        // Some known values
        int N = Lattice*Lattice; // Number of states/elements
        double E_empty = 0; // Initial energy
        double M_empty = 0; // Initial magnetism
        double J = 1.0; // Coupling constant = 1
        double beta =1/(T*k_b);
        // Creating an 'empty' matrix (filled with zeros)
        arma::mat S_empty(Lattice, Lattice);
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
        // Parallel loop over T
       
       

        /*====================================*/
        /*~~~~  Parallelized loop over T  ~~~~*/
        /*====================================*/
        #pragma omp for
        for(int t = 0; t < n_T; t++){
            // 
            // Inner loop with some computation for the given T
            //
            
            // beta avhenger av temp!!! Put verything that depends on 
            // temperature and beta in this very loop!!!

            // MCMC-loop
            for (int no = 0; no < n_cycles; no++)
            {
            // Cycle number
            int cycle = no + 1;
            // Empty vector for the new_E below
            vector<double> list_Es_new = {};

            arma::mat S_MC      = MCMC_s.random_spinnergal(S,T,L,N,E,M,beta,plusone,minusone);
            //cout << "\nEnergy before tot_energyboi" << E;
        

            // Expected energy per spin calculations       
            double new_E 	    = MCMC_s.tot_energyboi(S_MC, L, E, T,plusone,minusone,list_Es_new);
            sum_E               += new_E;   //energy pr lattice (for avg)
            sum_e     	        += new_E/N; //energy pr spin (for avg)
            double E_avg        = sum_E/cycle; // mean energy///////////////////////////////
            double e_avg 	    = sum_e/cycle; // mean energy per spin
            double e_avg_sqrd   = sum_e*sum_e/cycle; // mean energy squared
        //    cout << "\ne_avg: " << E_avg << " exp_E: " << exp_E;
        
            // Expected magnetization pr spin calculations
            double new_M		= MCMC_s.tot_magnetboi_abs(S_MC, T, L, M);
            double exp_m	    = abs(new_M/N); // magnetization per spin
            // standard deviation 
            double std_dev_e = e_avg_sqrd - e_avg*e_avg; 

            // Writing the values in the datafile
            datafile << cycle << " "<< e_avg << " " << exp_m << endl;
            }// end of for loop(no = cycles)
        }// end of parallel loop

//        datafile.close(); //file is finished
    }// end of parallel region

    // end of main function
    return 0;
}
