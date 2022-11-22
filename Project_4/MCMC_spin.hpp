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

// Returns a random (int) number between (and including) start and stop
int MC_generator(int start, int stop);

// std::vector<int> plus_minus_bois(int L); //-> tuple<vector<int>, vector<int>>;

// Spins the elements in lattice S
arma::mat spinnerboi(arma::mat S, int L);

// Calculates the total energy of the lattice
double tot_energyboi(arma::mat S, int L, double E, double T, std::vector<int> plusone, std::vector<int> minusone, std::vector<double>& list_Es);

// Creates a vector with energy for each atom
//std::vector<double> energy_listboi(arma::mat S, int L, std::vector<int> plusone, std::vector<int> minusone);

// Calculates the total magnetism of the lattice
double tot_magnetboi(arma::mat S, double T, int L, double M);
double tot_magnetboi_abs(arma::mat S, double T, int L, double M_abs);

// Finding the Boltzman
double boltzman_factors(double beta,int delta_E);

// Spins a single random atom in the lattice
arma::mat random_spinnergal(arma::mat& S, double T, int L, int N, double& E, double& M, double beta, std::vector<int> plusone, std::vector<int> minusone);


}; // Classes always end with };

#endif
