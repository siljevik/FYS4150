#ifndef __MCMC_spin_hpp__   
#define __MCMC_spin_hpp__

class MCMC_spin 
{
private:
  // Declaration of variables only accessible from within the class

public:
  // Declaration of constructors, e.g.

// variables

double T; // temperature
double TC; // critical temperature

double kb; // Boltzmann's constant
arma::mat S; // LxL matrix with spin states

/*
  MCMC_spin(double T_in) 
  {
    T = T_in;
  }
*/
///// Functions
  // For Problem 4 we want to declear the following functions:

  // memberfunctions for generating S
  void rand_spin_matrix();

   // Initialize the omp_rng_container according to the number of threads
  void initialize_omp_rng_container(unsigned int base_seed=-1);
  void get_random_int_m1_p1();


// Basic but highly necessary functions for later use
  void energy(amra::mat S); // energy dependant on lattice L
  void partition_Z(double beta, double J); // returns partition Z
  void prob_func(double Z, double beta, ); // probability function



  // define the expactationvalues of E and M
  void exval_E(); // return the expected energy
  void exval_E_sqr();

  void exval_B(); // return expected magnetization in absolute value(else it will be 0)
  void exval_B_sqrd();
  // expected energy (<epsilon>)
  void epsilon(int N, arma::vec E); //return expected energy pr spin 
  // expected magnetization (<|m|>)
  void abs_mag(int N, arma::vec);// return expected magnetization pr spin


  // Sepcific heat capacity (C_V)
  void spec_heatcap(int N, double T, double k_b/* <E^2>, <E>*/);
  // Susceptibility (chi). Normalized to number of spins
  void norm_sus(int N, double T, double k_b /*<M^2>, <M>*/);

  MCMC_spin(arguments);

  // Declaration of destructors, copy constructors, ...

  // Declarations of other class methods, e.g.
  void some_function(arguments);

}; // <-- Note that class bodies end with a semicolon!

#endif