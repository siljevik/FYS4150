#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <chrono>

/////////////////////////////
//#pragma                      WHAT IS THIS
/////////////////////////////
#ifndef __analytical_hpp__
#define __analytical_hpp__


class analytical 
{
private:
// Nothing to add?

public:

// Partition function
double part_func(double J, double beta);
/*
// Expected total energy
double exp_tot_E(double J, double beta, double Z);
double exp_tot_E_sqrd(double J, double beta, double Z);

// Expected total magentization
double exp_tot_M(double J, double beta, double Z);
double exp_tot_M_sqrd(double J, double beta, double Z);

//Specific heat capacity
double spec_heat_cap(int N, double J, double beta, double kb, double T, double exp_E, double exp_EE);
//Susceptibility
double sus_chi(int N, double J, double beta, double kb, double T, double exp_M, double exp_MM);
*/
}; // Classes always end with };

#endif