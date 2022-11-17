#ifndef __analytical_hpp__
#define __analytical_hpp__


class analytical{
    
    public:
    // Dette gjør rare ting med koden i analytical.cpp. Why???????
    /*int L;
    int N;
    double const kb;
    double T;
    double J;*/
    //double beta; try implementing this in each function?

    // Partition function
    double part_func(double J, double beta);

    // Expected total energy
    double exp_tot_E(double J, double beta, double Z);
    double exp_tot_E_sqrd(double J, double beta, double Z);

    // Expected total magnetism
    double exp_tot_M(double J, double beta, double Z);
    double exp_tot_M_sqrd(double J, double beta, double Z);

    //Specific heat capacity
    double spec_heat_cap(double N, double kb, double T, double exp_E, double exp_EE);
    //Susceptibility
    double sus_chi(double N, double kb, double T, double exp_M, double exp_MM);

}; // Classes always end with };

#endif