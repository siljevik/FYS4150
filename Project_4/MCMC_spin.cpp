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




// Analytical expression of Z for 2D lattice (2x2)
//ouble MCMC_spin::partition_Z(double beta, double J){
  //for problem 4 we have the partition function
  //  double Z = 4 * cosh(8*beta*J) + 12;
  //  return Z;
//}

