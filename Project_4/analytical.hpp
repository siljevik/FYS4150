#ifndef __analytical_hpp__
#define __analytical_hpp_

double const kb;
double T;
double J;
double beta;

// Partition function
double part_func(double J, double beta);

// Expected total energy
double exp_tot_E(double J, double beta);
double exp_tot_E_sqrd(double J, double beta);

// Expected total magentization
double exp_tot_M(double J, double beta);
double exp_tot_M_sqrd(double J, double beta);

//Specific heat capacity
double spec_heat_cap(int N, double J, double beta, double kb, double T);
//SUSceptibility
double sus_chi(int N, double J, double beta, double kb, double T);

#endif