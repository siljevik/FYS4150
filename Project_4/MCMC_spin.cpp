#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <random>
#include <vector>
#include <map>


#include "MCMC_spin.hpp"
#include "analytical.hpp"

using namespace std;



/*===========================================*/
/*~~~~~~     Monte Carlo Generator     ~~~~~~*/
/*===========================================*/
// Creates a random number between (and including) start and stop
int MCMC_spin::MC_generator(int start, int stop)
{
  // The random number generator from example:
  // https://github.com/anderkve/FYS3150/blob/master/code_examples/random_number_generation/main_basics.cpp
  unsigned int seed = chrono::system_clock::now().time_since_epoch().count();
  mt19937 generator;
  generator.seed(seed);
  uniform_int_distribution<int> my_01_pdf(start,stop);
  // Picking random value between (and including) start and stop
  int random = my_01_pdf(generator);
  return random;
}



/*===========================================*/
/*~~~~~~ Generator of spins in lattice ~~~~~~*/
/*===========================================*/
// This is only for filling up the empty matrix of size 
// LxL with random spins (the first time only)
arma::mat MCMC_spin::spinnerboi(arma::mat S, int L)
{
  // The spins are either up (+1) or down (-1), and constant 
  // since that won't change. So this is a vector we will use 
  // for the random generating of spin up or down.
  const vector<int> spinlist = {-1,+1};
  // For each element in matrix S
  for(int i = 0; i<L; i++){
    for(int j = 0; j<L; j++){
      // Replacing the value in the matrix with the random spin
      int rand_spin = spinlist[MCMC_spin::MC_generator(0,1)];
      S(i,j) = rand_spin;
    }
  }
  return arma::mat (S);
}



/*==========================================*/
/*~~~~~ Total energy of the 2D lattice ~~~~~*/
/*==========================================*/
double MCMC_spin::tot_energyboi(arma::mat S, int L, double E, double T,vector<int> plusone, vector<int> minusone, vector<double>& list_Es)
{
  // Looping though the Lattice
  for(int i = 0; i<L; i++)
  { 
    // Making the plusone and minusone in the x-directions to 
    // speed up the code a little bit
    int poi = plusone[i];
    int moi = minusone[i];
    
    for(int j = 0; j<L; j++){
      // Placement of the state we are doing calculations for
      int E_ij = S(i,j);

      // Surrounding atoms
      int E_o = S(i,minusone[j]);   //over
      int E_u = S(i,plusone[j]);   //under
      int E_v = S(moi,j);   //v for left (NO)
      int E_h = S(poi,j);    //h for right (NO)
      
      // Adding all energies to the total energy
      E += (-E_ij*(E_o + E_u + E_v + E_h))/2; // divided by 2 to correct for the doublecounting
      // Adding the calculated energy into the list_Es
      list_Es.push_back(E);
    }
  }
  //cout << "Energy list: " << tot_energy_pr_atom_list;
  return E;
}



/*=============================================*/
/*~~~~~ Total magnetism of the 2D lattice ~~~~~*/
/*=============================================*/
double MCMC_spin::tot_magnetboi(arma::mat S, double T, int L, double M)
{
  // Looping though the lattice
  for(int i = 0; i<L; i++)
  { 
    for(int j = 0; j<L; j++)
    {
      // Each spin
      int M_ij = S(i,j);
      // Adding all spins to the total magnetism
      M += M_ij;
    }
  }
  return M;
}


double MCMC_spin::tot_magnetboi_abs(arma::mat S, double T, int L, double M_abs)
{
  // Looping though the lattice
  for(int i = 0; i<L; i++)
  { 
    for(int j = 0; j<L; j++)
    {
      // Each spin
      int M_ij = abs( S(i,j));
      // Adding all spins to the total magnetism
      M_abs += abs(M_ij);
    }
  }
  return M_abs;
}

double MCMC_spin::boltzman_factors(double beta,int delta_E)
{
  // Make it into map????
  map<int,double> boltzman_;
  boltzman_[8] = exp(-beta*8); // our p_sT
  boltzman_[4] = exp(-beta*4);
  boltzman_[0] = exp(-beta*0);
  boltzman_[-4] = exp(-beta*(-4));
  boltzman_[-8] = exp(-beta*(-8));
  double boltzman = boltzman_[delta_E];
  return boltzman;
}



/*==================================================*/
/*~~~~~~ Generator of random spins in lattice ~~~~~~*/
/*==================================================*/    // & forran variabel oppdaterer også i main()

// Implementing a Metropolis algorithm as random_spinnergal
arma::mat MCMC_spin::random_spinnergal(arma::mat& S, double T, int L, int N, double& E, double& M, double beta, vector<int> plusone, vector<int> minusone)
{
  // Getting random indices for matrix S (= S(x,y))
  int x = MCMC_spin::MC_generator(0,(L-1));
  int y = MCMC_spin::MC_generator(0,(L-1));
  // Random 0 or 1
  int r = MCMC_spin::MC_generator(0,1);
  
  // Sum of the surrounding particles
  int surr_sum = S(minusone[x],y)+S(plusone[x],y)+S(x,minusone[y])+S(x,plusone[y]); 
  // Difference between initial and final
  int delta_E = 2*(-S(x,y))*surr_sum; 

  double boltzman = MCMC_spin::boltzman_factors(beta,delta_E);

//  cout <<"\n delta_E: "<<delta_E<<" boltzman: "<<boltzman<<"\n";
  
  // Should the spin be flipped? (Mac: alt+7 = |, also, here || means or)
  if (delta_E <= 0 || r <= boltzman)
  {
    // Flipping the spin here
    S(x,y) = - S(x,y);
    // Updating energy and magnetism
    // E = MCMC_spin::tot_energyboi(S,L,E,T);
    E += (double) delta_E;
    M += (double) 2*S(x,y);
  }
  return arma::mat (S);
}



////////////////////////////////////////////
// NOKE ME FEKK AV EIN GRUPPELÆRAR
//prob_after/prob_initial
//if A is not r 0 or 1 (generated randomly)
////////////////////////////////////////////



double MCMC_spin::prob_func(double beta, double E_before, double E_after, double Z)
{
  // want to calculate the Boltzmann prob. dist.
  double p_sT = (1/Z)*exp(-beta*(E_after - E_before) );
  return p_sT;
}


