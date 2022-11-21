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

/*
vector<int> MCMC_spin::plus_minus_bois(int L)
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
  return minusone, plusone;
}
*/

/*===========================================*/
/*~~~~~~ Generator of spins in lattice ~~~~~~*/
/*===========================================*/
// Since we want to return the updated vector

// This is only for making the first matrix
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
double MCMC_spin::tot_energyboi(arma::mat S, int L, double E, double T)
{
  // Creating empty vectors to look up indexes so we can look
  // at neighbours of the states (in case of 'border'-states)
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
      // Placement of the state we are doing calculations for
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

      // push the calculated energy into the tot_energy list
      tot_energy_pr_atom_list.push_back(E);
    }
  }
  //cout << "Energy list: " << tot_energy_pr_atom_list;
  return tot_energy_pr_atom_list;
}



/*=============================================*/
/*~~~~~ Total magnetism of the 2D lattice ~~~~~*/
/*=============================================*/
double MCMC_spin::tot_magnetboi(arma::mat S, double T, int L, double M)
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
arma::mat MCMC_spin::random_spinnergal(arma::mat& S, double T, int L, int N, double& E, double& M, double beta)
{
    // The random number generator from example:
    // https://github.com/anderkve/FYS3150/blob/master/code_examples/random_number_generation/main_basics.cpp
    unsigned int seed = chrono::system_clock::now().time_since_epoch().count();
    mt19937 generator;
    generator.seed(seed);
    uniform_int_distribution<int> my_01_pdf(0,(L-1) );
    uniform_int_distribution<int> my_02_pdf(0,(L-1) ); // L-1 to keep within range
    uniform_int_distribution<int> my_03_pdf(0,1);

    // Picking random indicies
    int x = my_01_pdf(generator);
    int y = my_02_pdf(generator);

    // Random 0 or 1
    int r = my_03_pdf(generator);

    /*============================================*/
    /*= Repeating the success from tot_energyboi =*/
    /*============================================*/

    // Creating empty vectors to look up indexes so we can look
    // at neighbours of the states (in case of 'border'-atoms)
    vector<int> plusone{};
    vector<int> minusone{};

    // Then fill up plusone (e.g. [1,2,3, ... , L-1, 0])
    for (int i = 0; i<L-1; i++){
      plusone.push_back(i+1);
    }// end of for loop(i)

    plusone.push_back(0);

    // Fill up minusone (e.g. [L-1,0,1,2, ... , L-2])
    minusone.push_back(L-1);
    for (int i = 0; i<L-1; i++){
      minusone.push_back(i);
    }// end of for loop(i)

    // Change of energy
    for(int j = 0; j < L; j++)
    {
    // Same as earlier, we speed up the code by looking in the "x-direction"
    // of the plu- and minusone vectors
    int poi = plusone[j];
    int moi = minusone[j];

    for(int k = 0; k < L; k++)
    {

    int surr_sum = S(moi,y)+S(poi,y)+S(x,minusone[k])+S(x,plusone[k]); // Sum of the surrounding particles
    int delta_E = 2*S(x,y)*surr_sum; // Difference between initial and final
    double boltzman = MCMC_spin::boltzman_factors(beta,delta_E);

// Should the spin be flipped? (Mac: alt+7 = |, also, here || means or)
    if (delta_E <= 0 || r <= boltzman)
    {
      // Flipping the spin here
      S(x,y) = - S(x,y);
      // Updating energy and magnetism
      //E = MCMC_spin::tot_energyboi(S,L,E,T);
	E = delta_E;
      M = MCMC_spin::tot_magnetboi(S,L,M,T);

    }// end of if statement


    } // end of for loop(k)
  }// end of for loop(j)

  return arma::mat (S);

}// end of spinnergal
 // expectation values / number of MC cycles




//prob_after/prob_initial
//if A is not r 0 or 1 (generated randomly)





double MCMC_spin::prob_func(double beta, double E_before, double E_after, double Z)
{
  // want to calculate the Boltzmann prob. dist.
  double p_sT = (1/Z)*exp(-beta*(E_after - E_before) );
  return p_sT;
}


