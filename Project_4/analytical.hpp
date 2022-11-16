#ifndef __analytical_hpp__
#define __analytical_hpp__


class Analytical{

    public:
    // Dette gjør rare ting med koden i analytical.cpp. Why???????
    int L;
    int N;
    double const kb;
    double T;
    double J;
    //double beta; try implementing this in each function?

<<<<<<< HEAD
// Partition function
double part_func(double J, double beta);

// Expected total energy
double exp_tot_E(double J, double beta, double Z);
double exp_tot_E_sqrd(double J, double beta, double Z);
=======
>>>>>>> 2f0eed05e0d2bcdd0814089a51d51c08d51d8de9

    //constructor
    analytical(double J_in, double kb_in, double T_in);

<<<<<<< HEAD
//Specific heat capacity
double spec_heat_cap(int N, double J, double beta, double kb, double T, double exp_E, double exp_EE);
//Susceptibility
double sus_chi(int N, double J, double beta, double kb, double T, double exp_M, double exp_MM);

}; // Classes always end with };
=======
    // Partition function
    double part_func(double J_in, double kb_in, double T_in);
>>>>>>> 2f0eed05e0d2bcdd0814089a51d51c08d51d8de9

    // Expected total energy
    double exp_tot_E(double J_in, double kb_in, double T_in);
    double exp_tot_E_sqrd(double J_in, double kb_in, double T_in);

    // Expected total magentization
    double exp_tot_M(double J_in, double kb_in, double T_in);
    double exp_tot_M_sqrd(double J_in, double kb_in, double T_in);

    //Specific heat capacity
    double spec_heat_cap(int N_in, double J_in, double kb_in, double T_in);
    //SUSceptibility
    double sus_chi(int N_in, double J_in, double kb_in, double T_in);

};
#endif