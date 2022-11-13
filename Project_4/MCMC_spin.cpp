#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <random>

#include "MCMC_spin.hpp"

using namespace std;


MCMC_spin::MCMC_spin()
{
  // ...
};

void MCMC_spin::spinnerboi(arma::mat S)
{
  // The spins are either up (+1) or down (-1), 
  // and constant since that won't change.
  const vector<int> spinlist(-1,+1);

  // For each element in matrix S
  for(int i = 0; i<L-1; i++)
  {
    for(int j = 0; j<L-1; j++)
    {
      // Choosing a random spin
      int a = rand() % 100;
      // Replacing the value in the matrix with the random spin
      int S(i,j) = a;
    }
  }
}



// Totl energy of the 2D lattice 
double MCMC_spin::tot_energy(arma::mat S, int L)
{
  // Creating empty vectors to look up indexes so we can look
  // at neighbours of the states (in case of 'border'-atoms)
  vector<int> plusone;
  vector<int> minusone;
  
  // Then fill up plusone (e.g. [1,2,3, ... , L-1, 0])
  for (int i = 0; i<L-1; i++){
    plusone.pushback(i+1);
  }
  plusone.pushback(0);
  
  // Fill up minusone (e.g. [L-1,0,1,2, ... , L-2])
  minusone.pushback(L-1);
  for (int i = 0; i<L-1; i++){
    minusone.pushback(i);
  }
  
  // Creating an E integer for the total energy of the system
  int E = 0; 
  
  // Looping though the 
  for(int i = 0; i<L-1; i++)
  { 
    // Making the plusone and minusone in the x-directions to 
    // speed up the code a little bit
    int poi = plusone[i];
    int moi = minusone[i];
    
    for(int j = 0; j<L-1; j++){
      // Placement of the atom we are doing calculations for
      int E_ij = S(i,j);
      
      // Surrounding atoms
      int E_o = S(i,minusone[j]);   //over
      int E_u = S(i,plusone[j]);   //under
      int E_v = S(moi,j);   //v for left (NO)
      int E_h = S(poi,j)    //h for right (NO)
      
      // Adding all energies to the total energy
      int E += -E_ij*(E_o + E_u + E_v + E_h);
      // return ze sum of spin divided by 2 to correct for the doublecounting
      return E/2;
    }
  }
}






// Analytical expression of Z for 2D lattice (2x2)
//ouble MCMC_spin::partition_Z(double beta, double J){
  //for problem 4 we have the partition function
  //  double Z = 4 * cosh(8*beta*J) + 12;
  //  return Z;
//}

