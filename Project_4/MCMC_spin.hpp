#ifndef __MCMC_spin_hpp__   
#define __MCMC_spin_hpp__

class MCMC_spin 
{
private:
// Nothing to add?

public:

// Spins the elements in lattice S
void spinnerboi(arma::mat S);

// Calculates the total energy of the lattice
void tot_energy(arma::mat S, int L);


}; // Classes always end with };

#endif