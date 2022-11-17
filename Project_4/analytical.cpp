// Here we implement all the analytical expressions for later use
#include <cmath>
#include <random>
#include <vector>
#include <math.h>


#include "analytical.hpp"
//#include "MCMC_spin.hpp"

using namespace std;

/*~~~~~~~ Implemented analytical expressions ~~~~~~~*/

// fix the damn variables!!!

/*==================================*/
/*~~~~~~ Partition function Z ~~~~~~*/
/*==================================*/

// Analytical expression og the partition function, Z

double analytical::part_func(double J, double beta)
{
    double c = 8*beta*J;
    // Expression Z
    double Z = 4*cosh(c) + 12;
    return Z;
}


/*==================================*/
/*~~~~~ Expected energy values ~~~~~*/
/*==================================*/

// Expected total energy
double analytical::exp_tot_E(double J, double beta,double Z)
{
    double c = 8*beta*J;
    double invZ = 1/Z;
    // Expected energy, <E>
    double exp_E = (-J)*8*sinh(c)*invZ;
    return exp_E;
}

// Expected squared total energy, <E^2>

double analytical::exp_tot_E_sqrd(double J, double beta, double Z)
{
    double c = 8*beta*J;
    double JJ = J*J;
    double invZ = 1/Z;
    // <E^2>
    double exp_EE = JJ*64*sinh(c)*invZ;
    return exp_EE; 
}


/*==================================*/
/*~~~~~ Expected magnetization ~~~~~*/
/*==================================*/
// Expected total magentization
double analytical::exp_tot_M(double J, double beta, double Z)
{
    double c = 8*beta*J;
    double invZ = 1/Z;
    double top = 8*exp(c) + 6;
    // <|M|>
    double exp_M = top*invZ;
    return exp_M;
}

// Expected squared magnetization

double analytical::exp_tot_M_sqrd(double J, double beta, double Z)
{
    double c = 8*beta*J;
    double invZ = 1/Z;
    double top = 32*exp(c) + 8;
    // <M^2>
    double exp_MM = top*invZ;
    return exp_MM;
}


/*==================================*/
/*~~~~~ Specific heat capacity ~~~~~*/
/*==================================*/
//Specific heat capacity, CV, normalized to number of spins, N:
double analytical::spec_heat_cap(int N, double J, double beta, double kb, double T, double exp_E, double exp_EE)
{
    double broek = ( 1/N )*( 1/(kb*T) );
    double variance_E = exp_EE - (exp_E*exp_E);
    // CV
    double CV = broek*variance_E;
    return CV;
}


/*==================================*/
/*~~~~~~~~~ Susceptibility ~~~~~~~~~*/
/*==================================*/
//Susceptibility, chi, normailzed to number of spins, N:
double analytical::sus_chi(int N, double J, double beta, double kb, double T, double exp_M, double exp_MM)
{
    double broek = ( 1/N )*( 1/(kb*T) );
    double variance_M = exp_MM - (exp_M*exp_M);
    // chi
    double chi = broek*variance_M;
    return chi;
}
