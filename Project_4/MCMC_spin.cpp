#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <random>
#include <vector>


#include "MCMC_spin.hpp"
#include "analytical.hpp"

using namespace std;


/*===========================================*/
/*~~~~~~ Generator of spins in lattice ~~~~~~*/
/*===========================================*/
// Since we want to return the updated vector
arma::mat MCMC_spin::spinnerboi(arma::mat S, int L)
{
  // The spins are either up (+1) or down (-1), and constant 
  // since that won't change. So this is a vector we will use 
  // for the random generating of spin up or down.
  const vector<int> spinlist = {-1,+1};

  // The random number generator from example: 
  // https://github.com/anderkve/FYS3150/blob/master/code_examples/random_number_generation/main_basics.cpp 
  unsigned int seed = chrono::system_clock::now().time_since_epoch().count();
  mt19937 generator;
  generator.seed(seed);
  uniform_int_distribution<int> my_01_pdf(0,1);

  // For each element in matrix S
  for(int i = 0; i<L; i++)
  {
    for(int j = 0; j<L; j++)
    {
      // Replacing the value in the matrix with the random spin
      S(i,j) = spinlist[my_01_pdf(generator)];
    }
  }
  return arma::mat (S);
}


/*==========================================*/
/*~~~~~ Total energy of the 2D lattice ~~~~~*/
/*==========================================*/
double MCMC_spin::tot_energyboi(arma::mat S, int L, int E)
{
  // Creating empty vectors to look up indexes so we can look
  // at neighbours of the states (in case of 'border'-atoms)
  vector<int> plusone{};
  vector<int> minusone{};
  
  // Then fill up plusone (e.g. [1,2,3, ... , L-1, 0])
  for (int i = 0; i<L-1; i++){
    plusone.push_back(i+1);
  }
  plusone.push_back(0);
  // Fill up minusone (e.g. [L-1,0,1,2, ... , L-2])
  minusone.push_back(L-1);
  for (int i = 0; i<L-1; i++){
    minusone.push_back(i);
  }
  
  // Looping though the Lattice
  for(int i = 0; i<L; i++)
  { 
    // Making the plusone and minusone in the x-directions to 
    // speed up the code a little bit
    int poi = plusone[i];
    int moi = minusone[i];
    
    for(int j = 0; j<L; j++){
      // Placement of the atom we are doing calculations for
      int E_ij = S(i,j);
      
      // Surrounding atoms
      int E_o = S(i,minusone[j]);   //over
      int E_u = S(i,plusone[j]);   //under
      int E_v = S(moi,j);   //v for left (NO)
      int E_h = S(poi,j);    //h for right (NO)
      
      // Adding all energies to the total energy
      E += (-E_ij*(E_o + E_u + E_v + E_h))/2; // divided by 2 to correct for the doublecounting
    }
  }
  //cout << "Energy list: " << tot_energy_pr_atom_list;
  return E;
}

/*==========================================================*/
/*~~~~~ List of energy for each atom in the 2D lattice ~~~~~*/
/*==========================================================*/
vector<double> MCMC_spin::energy_listboi(arma::mat S, int L)
{
  // Creating empty vectors to look up indexes so we can look
  // at neighbours of the states (in case of 'border'-atoms)
  vector<int> plusone{};
  vector<int> minusone{};

  // For the degeneracy
  vector<double> tot_energy_pr_atom_list{};

  // Energy per atom
  double E;
  
  // Then fill up plusone (e.g. [1,2,3, ... , L-1, 0])
  for (int i = 0; i<L-1; i++){
    plusone.push_back(i+1);
  }
  plusone.push_back(0);
  // Fill up minusone (e.g. [L-1,0,1,2, ... , L-2])
  minusone.push_back(L-1);
  for (int i = 0; i<L-1; i++){
    minusone.push_back(i);
  }
  
  // Looping though the Lattice
  for(int i = 0; i<L; i++)
  { 
    // Making the plusone and minusone in the x-directions to 
    // speed up the code a little bit
    int poi = plusone[i];
    int moi = minusone[i];
    
    for(int j = 0; j<L; j++){
      // Placement of the atom we are doing calculations for
      int E_ij = S(i,j);
      
      // Surrounding atoms
      int E_o = S(i,minusone[j]);   //over
      int E_u = S(i,plusone[j]);   //under
      int E_v = S(moi,j);   //v for left (NO)
      int E_h = S(poi,j);    //h for right (NO)
      
      // Adding all energies to the total energy
      E = (-E_ij*(E_o + E_u + E_v + E_h))/2; // divided by 2 to correct for the doublecounting
      tot_energy_pr_atom_list.push_back(E);
    }
  }
  //cout << "Energy list: " << tot_energy_pr_atom_list;
  return tot_energy_pr_atom_list;
}


/*=============================================*/
/*~~~~~ Total magnetism of the 2D lattice ~~~~~*/
/*=============================================*/
double MCMC_spin::tot_magnetboi(arma::mat S, int L, int M)
{
  // Looping though the lattice
  for(int i = 0; i<L; i++)
  { 
    for(int j = 0; j<L; j++){
      // Each spin
      int M_ij = S(i,j);
      // Adding all spins to the total magnetism
      M += M_ij;
    }
  }
  return M;
}

arma::mat MCMC_spin::single_spinnergal(arma::mat S, int L, double boltzman_value)
{
  for (int i = 0; i < N; i++)
  {
    // The random number generator from example: 
    // https://github.com/anderkve/FYS3150/blob/master/code_examples/random_number_generation/main_basics.cpp 
    unsigned int seed = chrono::system_clock::now().time_since_epoch().count();
    mt19937 generator;
    generator.seed(seed);
    uniform_int_distribution<int> my_01_pdf(0,L-1);
    uniform_int_distribution<int> my_02_pdf(0,L-1); // L-1 to keep within range
    uniform_int_distribution<double> my_03_pdf(0,1);
    // Picking random indicies
    int x = my_01_pdf(generator);
    int y = my_02_pdf(generator);
    double r = my_03_pdf(generator);
    
    // Change of energy
    double surr_sum = S(x-1,y)+S(x+1,y)+S(x,y-1)+S(x,y+1);
    double delta_E = 2*S(x,y)*surr_sum; // Difference between initial and final
    double boltzman_n = delta_E

    // Should the spin be flipped? (Mac: alt+7 = |, also, here || means or)
    if (delta_E < 0 || r <= MCMC_spin::boltzman_factors(boltzman_n))
    {
      // Updating energy and magnetism
      E = MCMC_s.tot_energyboi(S2,L,E);
      M = MCMC_s.tot_magnetboi(S2,L,M);
      S(x,y) = - S(x,y); // Flipping the one random spin here
    }
  }
  return arma::mat (S);
}

double MCMC_spin::boltzman_factors(boltzman_n)
{
  // Do the map thing
  vector<double> boltzman_values = {exp(beta*8),exp(beta*4),exp(beta*0),exp(beta*(-4)),exp(beta*(-8))};
  double boltzman_value = boltzman_values[boltzman_n];
  
  boltzman_[8] = exp(beta*8);
  boltzman_[4] = exp(beta*4);
  boltzman_[0] = exp(beta*0);
  boltzman_[-4] = exp(beta*(-4));
  boltzman_[-8] = exp(beta*(-8));
  return boltzman_value

}

prob_after/prob_initial
if A is not r 0 or 1 (generated randombly)


double MCMC_spin::prob_func(double beta, double E_before, double E_after, double Z)
{
  // want to calculate the Boltzmann prob. dist.
  double p_sT = (1/Z)*exp(-beta*(E_after - E_before) );
  return p_sT;
}

/*
// Inspired by the code at page 150 from the course book/pdf
void MC_sampling(double Z, double beta, double E_s)
//int initial_n_particles, int max_time,int number_cycles, double decay_probability,int *ncumulative)
{
  // Given a system temperature , the probability for the system state  is given by the Boltzmann distribution:
  double p_sT = (1/Z)*exp(-(beta)*E_s); // E_s = E2 in main, total energy

  // Probability distribution (prob_distr)
  double expected_E = (E_s)*p_sT; 


    int cycles, time, np, n_unstable, particle_limit;
    long idum;
    idum=-1; // initialise random number generator
    // loop over monte carlo cycles
    // One monte carlo loop is one sample
    for (cycles = 1; cycles <= number_cycles; cycles++)
    {
        n_unstable = initial_n_particles;
        // accumulate the number of particles per time step per trial
        ncumulative[0] += initial_n_particles;
        
        // loop over each time step
        for (time=1; time <= max_time; time++)
        {
            // for each time step, we check each particle
            particle_limit = n_unstable;
            for ( np = 1; np <= particle_limit; np++) 
            {
                if( ran0(&idum) <= decay_probability) 
                {
                    n_unstable=n_unstable-1;
                }
            }
            ncumulative[time] += n_unstable;
        }
    }
}


*/

