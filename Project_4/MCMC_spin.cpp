#include "MCMC_spin.hpp"
#include "rand_state.cpp"

// Definitions of constructors
/*
MCMC_spin::MCMC_spin(arguments)
{
  // ...
}
*/


 // need to fill the S matrix with random values -1 or +1
  // use get_random_int_0_1 as base to find random number -1 or +1
  /*
  - use a similiar function to fill an array with random numbers -1 or +1
  - fill the S matrix with the values
  */

// call get_random_int_m1_p1() which returns random int r


arma::mat rand_spin_matrix(int L, int r){ 
  arma::mat S(L,L,fill::r); // define the matrix S
    
 
  return S;
};

// Totl energy of the 2D lattice 
double MCMC_spin::energy(arma::mat S, int L){
// Implement the sum over all spin pairs <kl>

/* for all pairs <kl> sum sk*sl N times and divide by two to 
 correct for the doublecounting */

// Must create empty vectors to look up indexilicious neighbours of the states
arma::vec plusone = [];
arma::vec minusone = [];

// Then fill up plusone
for (int i = 0; i<L-1; i++);{
  plusone.pushback(i+1);
};
plusone.pushback(0);

// Fill up minusone
minusone.pushback(L-1)
for (int i = 0; i<L-1; i++);{
  minusone.pushback(i);
};

// use the indexing to check the neighbours of the spin-state.
int E = 0;
// run through the mf indexes 
// sum sk*sl; Bruk modulus for å sjekke
for(int i = 0; i<L-1; i++);{ //x-direction
    int poi = plusone[i];
    int moi = minusone[i];

    for(int j = 0; j<L-1; j++){
      // Plassering, atomet vi ser på
      int E_ij = S(i,j);

      // Surrounding atoms
      int E_o = S(i,minusone[j]);   //over
      int E_u = S(i,plusone[j]);   //under
      int E_v = S(moi,j);   //v for left (NO)
      int E_h = S(poi,j)    //h for right (NO)
      // Adding all energies
      E += -E_ij*(E_o + E_u + E_v + E_h);
    // return ze sum of spin divided by 2
    return E/2;
    };
};
}
// Analytical expression of Z for 2D lattice (2x2)
double MCMC_spin::partition_Z(double beta, double J){
  //for problem 4 we have the partition function
    double Z = 4 * cosh(8*beta*J) + 12;
    return Z;
}

// Definitions of copy constructor, ...

// Definitions of other class methods, e.g.
/* 
void MCMC_spin::some_function(arguments)
{
  // ...
}

*/