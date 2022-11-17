#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <random>

#ifndef __MCMC_spin_hpp__   
#define __MCMC_spin_hpp__

class MCMC_spin 
{
private:
// Nothing to add?

public:

// Spins the elements in lattice S
arma::mat spinnerboi(arma::mat S, int L);

// Calculates the total energy of the lattice
double tot_energyboi(arma::mat S, int L, int E);

// Creates a vector with energy for each atom
std::vector<double> energy_listboi(arma::mat S, int L);

// Calculates the total magnetism of the lattice
double tot_magnetboi(arma::mat S, int L, int M);

// Spins a single random atom in the lattice
arma::mat single_spinnergal(arma::mat S, int L);

double prob_func(double beta, double E_before, double E_after, double Z);

}; // Classes always end with };

#endif