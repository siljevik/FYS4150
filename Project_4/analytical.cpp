// Here we implement all the analytical expressions for later use
#include <cmath>
#include <random>
<<<<<<< HEAD
#include <vector>
#include <math.h>


#include "analytical.hpp"
//#include "MCMC_spin.hpp"

using namespace std;

=======

#include "analytical.hpp"
>>>>>>> 2f0eed05e0d2bcdd0814089a51d51c08d51d8de9

/*~~~~~~~ Implemented analytical expressions ~~~~~~~*/

// fix the damn variables!!!

/*==================================*/
/*~~~~~~ Partition function Z ~~~~~~*/
/*==================================*/

// Analytical expression og the partition function, Z
<<<<<<< HEAD
double analytical::part_func(double J, double beta)
{
    double c = 8*beta*J;
=======
double part_func(double J_in, double kb_in, double T_in){
    double beta = 1/(kb_in*T_in);
    double c = 8*beta*J_in;
>>>>>>> 2f0eed05e0d2bcdd0814089a51d51c08d51d8de9

    // Expression Z
    double Z = 4*cosh(c) + 12;

    return Z;
}


/*==================================*/
/*~~~~~ Expected energy values ~~~~~*/
/*==================================*/
<<<<<<< HEAD
// Expected total energy
double analytical::exp_tot_E(double J, double beta,double Z)
{
    double c = 8*beta*J;
=======

// Expected total energy
double exp_tot_E(double J_in, double kb_in, double T_in){
    double beta = 1/(kb_in*T_in);
    double c = 8*beta*J_in;
    double Z = 4*cosh(c) + 12;
>>>>>>> 2f0eed05e0d2bcdd0814089a51d51c08d51d8de9
    double invZ = 1/Z;

    // Expected energy, <E>
    double exp_E = (-J_in)*8*sinh(c)*invZ;
    return exp_E;
}

// Expected squared total energy, <E^2>
<<<<<<< HEAD
double analytical::exp_tot_E_sqrd(double J, double beta, double Z)
{
    double c = 8*beta*J;
    double JJ = J*J;
=======
double exp_tot_E_sqrd(double J_in, double kb_in, double T_in){
    double beta = 1/(kb_in*T_in);
    double c = 8*beta*J_in;
    double JJ = J_in*J_in;
    double Z = 4*cosh(c) + 12;
>>>>>>> 2f0eed05e0d2bcdd0814089a51d51c08d51d8de9
    double invZ = 1/Z;

    // <E^2>
    double exp_EE = JJ*64*sinh(c)*invZ;
    return exp_EE; 
}
<<<<<<< HEAD

=======
>>>>>>> 2f0eed05e0d2bcdd0814089a51d51c08d51d8de9

/*==================================*/
/*~~~~~ Expected magnetization ~~~~~*/
/*==================================*/

// Expected total magentization
<<<<<<< HEAD
double analytical::exp_tot_M(double J, double beta, double Z)
{
    double c = 8*beta*J;
=======
double exp_tot_M(double J_in, double kb_in, double T_in){
    double beta = 1/(kb_in*T_in);
    double c = 8*beta*J_in;
    double Z = 4*cosh(c) + 12;
>>>>>>> 2f0eed05e0d2bcdd0814089a51d51c08d51d8de9
    double invZ = 1/Z;
    double top = 8*exp(c) + 6;

    // <|M|>
    double exp_M = top*invZ;
    return exp_M;
}

// Expected squared magnetization
<<<<<<< HEAD
double analytical::exp_tot_M_sqrd(double J, double beta, double Z)
{
    double c = 8*beta*J;
=======
double exp_tot_M_sqrd(double J_in, double kb_in, double T_in){
    double beta = 1/(kb_in*T_in);
    double c = 8*beta*J_in;
    double Z = 4*cosh(c) + 12;
>>>>>>> 2f0eed05e0d2bcdd0814089a51d51c08d51d8de9
    double invZ = 1/Z;
    double top = 32*exp(c) + 8;

    // <M^2>
    double exp_MM = top*invZ;
    return exp_MM;
}
<<<<<<< HEAD

=======
>>>>>>> 2f0eed05e0d2bcdd0814089a51d51c08d51d8de9

/*==================================*/
/*~~~~~ Specific heat capacity ~~~~~*/
/*==================================*/
<<<<<<< HEAD
//Specific heat capacity, CV, normalized to number of spins, N:
double analytical::spec_heat_cap(int N, double J, double beta, double kb, double T, double exp_E, double exp_EE)
{
    double broek = ( 1/N )*( 1/(kb*T) );
=======

//Specific heat capacity, CV, normalized to number of spins, n:
double spec_heat_cap(int N_in, double J_in, double kb_in, double T_in){
    double broek = ( 1/N_in )*( 1/(kb_in*T_in) );
>>>>>>> 2f0eed05e0d2bcdd0814089a51d51c08d51d8de9
    double variance_E = exp_EE - (exp_E*exp_E);

    // CV
    double CV = broek*variance_E;
    return CV;
}
<<<<<<< HEAD


/*==================================*/
/*~~~~~~~~~ Susceptibility ~~~~~~~~~*/
/*==================================*/
//Susceptibility, chi, normailzed to number of spins, N:
double analytical::sus_chi(int N, double J, double beta, double kb, double T, double exp_M, double exp_MM)
{
    double broek = ( 1/N )*( 1/(kb*T) );
=======
/*==================================*/
/*~~~~~~~~~ Susceptibility ~~~~~~~~~*/
/*==================================*/

//SUSceptibility, chi, normailzed to number of spins, n:
double sus_chi(int N_in, double J_in, double kb_in, double T_in){
    double broek = ( 1/N_in )*( 1/(kb_in*T_in) );
>>>>>>> 2f0eed05e0d2bcdd0814089a51d51c08d51d8de9
    double variance_M = exp_MM - (exp_M*exp_M);

    // chi
    double chi = broek*variance_M;
    return chi;
<<<<<<< HEAD
}
=======
}
>>>>>>> 2f0eed05e0d2bcdd0814089a51d51c08d51d8de9
