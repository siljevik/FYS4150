// Here we implement all the analytical expressions for later use
#include <cmath>
#include <random>

#include "analytical.hpp"

/*~~~~~~~ Implemented analytical expressions ~~~~~~~*/

// fix the damn variables!!!

/*==================================*/
/*~~~~~~ Partition function Z ~~~~~~*/
/*==================================*/

// Analytical expression og the partition function, Z
double part_func(double J_in, double kb_in, double T_in){
    double beta = 1/(kb_in*T_in);
    double c = 8*beta*J_in;

    // Expression Z
    double Z = 4*cosh(c) + 12;

    return Z;
}


/*==================================*/
/*~~~~~ Expected energy values ~~~~~*/
/*==================================*/

// Expected total energy
double exp_tot_E(double J_in, double kb_in, double T_in){
    double beta = 1/(kb_in*T_in);
    double c = 8*beta*J_in;
    double Z = 4*cosh(c) + 12;
    double invZ = 1/Z;

    // Expected energy, <E>
    double exp_E = (-J_in)*8*sinh(c)*invZ;
    return exp_E;
}

// Expected squared total energy, <E^2>
double exp_tot_E_sqrd(double J_in, double kb_in, double T_in){
    double beta = 1/(kb_in*T_in);
    double c = 8*beta*J_in;
    double JJ = J_in*J_in;
    double Z = 4*cosh(c) + 12;
    double invZ = 1/Z;

    // <E^2>
    double exp_EE = JJ*64*sinh(c)*invZ;
    return exp_EE; 
}

/*==================================*/
/*~~~~~ Expected magnetization ~~~~~*/
/*==================================*/

// Expected total magentization
double exp_tot_M(double J_in, double kb_in, double T_in){
    double beta = 1/(kb_in*T_in);
    double c = 8*beta*J_in;
    double Z = 4*cosh(c) + 12;
    double invZ = 1/Z;
    double top = 8*exp(c) + 6;

    // <|M|>
    double exp_M = top*invZ;
    return exp_M;
}

// Expected squared magnetization
double exp_tot_M_sqrd(double J_in, double kb_in, double T_in){
    double beta = 1/(kb_in*T_in);
    double c = 8*beta*J_in;
    double Z = 4*cosh(c) + 12;
    double invZ = 1/Z;
    double top = 32*exp(c) + 8;

    // <M^2>
    double exp_MM = top*invZ;
    return exp_MM;
}

/*==================================*/
/*~~~~~ Specific heat capacity ~~~~~*/
/*==================================*/

//Specific heat capacity, CV, normalized to number of spins, n:
double spec_heat_cap(int N_in, double J_in, double kb_in, double T_in){
    double broek = ( 1/N_in )*( 1/(kb_in*T_in) );
    double variance_E = exp_EE - (exp_E*exp_E);

    // CV
    double CV = broek*variance_E;
    return CV;
}
/*==================================*/
/*~~~~~~~~~ Susceptibility ~~~~~~~~~*/
/*==================================*/

//SUSceptibility, chi, normailzed to number of spins, n:
double sus_chi(int N_in, double J_in, double kb_in, double T_in){
    double broek = ( 1/N_in )*( 1/(kb_in*T_in) );
    double variance_M = exp_MM - (exp_M*exp_M);

    // chi
    double chi = broek*variance_M;
    return chi;
}