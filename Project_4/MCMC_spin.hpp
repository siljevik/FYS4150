#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <random>
#include <map>

#ifndef __MCMC_spin_hpp__   
#define __MCMC_spin_hpp__

class MCMC_spin
{
private:
// Nothing to add?
// Variables
//double T;

public:
// Member variables - so the class know all the conditions
// double T; // Keeping track of the temperature


// std::vector<int> plus_minus_bois(int L); //-> tuple<vector<int>, vector<int>>;

// Spins the elements in lattice S
arma::mat spinnerboi(arma::mat S, int L);

// Calculates the total energy of the lattice
double tot_energyboi(arma::mat S, int L, double E, double T);

// Creates a vector with energy for each atom
std::vector<double> energy_listboi(arma::mat S, int L);

// Calculates the total magnetism of the lattice
double tot_magnetboi(arma::mat S, double T, int L, double M);

// Finding the Boltzman
double boltzman_factors(double beta,int delta_E);

// Spins a single random atom in the lattice
arma::mat random_spinnergal(arma::mat& S, double T, int L, int N, double& E, double& M, double beta);

// Finding the probability??????
double prob_func(double beta, double E_before, double E_after, double Z);


}; // Classes always end with };

#endif
