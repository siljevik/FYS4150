/////////////////
//
// Compile w/ : g++ double_test.cpp func/particle.cpp func/penningtrap.cpp func/analytical.cpp -std=c++11 -I include -o double_test.exe -larmadillo
//
// To execute file and save data (without interactions): ./double_test.exe > two_particle_data_no_int.txt
//
// To execute file and save data (with interactions): ./double_test.exe > two_particle_data_w_int.txt
//
/////////////////
#include <iostream>
#include <cmath>
#include <armadillo>
#include <iomanip> // needed for the setw() and setprecision()


// include classes
#include "header/particle.hpp"
#include "header/penningtrap.hpp"
#include "header/analytical.hpp"


using namespace std;

int main()
{
  int width = 18;
  int prec = 4;

  PenningTrap trap = PenningTrap(B0, V0, d);

  double q = 1.;	// charge of particle [e]
  double m = 40.078;  // mass of Ca+ ion [u]

  arma::vec r_1 = arma::vec{20.0, 0.0, 20.0}; // mu*m
  arma::vec r_2 = arma::vec{25.0, 25.0, 0.0}; // mu*m
  arma::vec v_1 = arma::vec{0.0, 25.0, 0.0}; // (mu*m) / (mu*s)
  arma::vec v_2 = arma::vec{0.0, 40.0, 5.0}; // (mu*m) / (mu*s)

  Particle p_1 = Particle(q, m, r_1, v_1);
  Particle p_2 = Particle(q, m, r_2, v_2);

  trap.add_particle(p_1);
  trap.add_particle(p_2);

  int t_tot = 50;   // total time [mu*s]
  double dt = pow(10,-3); // timestep dt
  int steps = t_tot / dt; // steps

  arma::vec t = arma::linspace(0, t_tot, steps);

  // Header naming the columns int the .txt-file
  cout << "Time"
       << "," << "x1"
       << "," << "x2"
       << "," << "y1"
       << "," << "y2"
       << "," << "z1"
       << "," << "z2"
       << "," << "vx1"
       << "," << "vx2"
       << "," << "vy1"
       << "," << "vy2"
       << "," << "vz1"
       << "," << "vz2"
       << endl;

  // First row in the .txt-file contains the initial values
  cout << t(0)
       << "," << r_1(0)
       << "," << r_2(0)
       << "," << r_1(1)
       << "," << r_2(1)
       << "," << r_1(2)
       << "," << r_2(2)
       << "," << v_1(0)
       << "," << v_2(0)
       << "," << v_1(1)
       << "," << v_2(1)
       << "," << v_1(2)
       << "," << v_2(2)
       << endl;

  for (int i = 1; i < steps; i++)
  {
	// Evolving the system w/ RK4 method and particle interaction
	// by changing between true or false we decide if there are particle interactions ##################### IMPORTANT DETAIL <3
	// true = particle interactions      --> for two_particle_data_w_int.txt
	// false = no particle interactions  --> for two_particle_data_no_int.txt
    	trap.evolve_RK4(dt, true);          // THIS IS WHERE WE CHANGE BETWEEN TRUE AND FALSE

    	double x1 = trap.particles[0].r[0];
	double x2 = trap.particles[1].r[0];
	double y1 = trap.particles[0].r[1];
    	double y2 = trap.particles[1].r[1];
    	double z1 = trap.particles[0].r[2];
    	double z2 = trap.particles[1].r[2];
    	double vx1 = trap.particles[0].v[0];
    	double vx2 = trap.particles[1].v[0];
    	double vy1 = trap.particles[0].v[1];
    	double vy2 = trap.particles[1].v[1];
    	double vz1 = trap.particles[0].v[2];
    	double vz2 = trap.particles[1].v[2];

    // Inserting the updated values for each row in the .txt-file
    cout << setprecision(prec) << t(i)
         << "," << setprecision(prec) << x1
         << "," << setprecision(prec) << x2
         << "," << setprecision(prec) << y1
         << "," << setprecision(prec) << y2
         << "," << setprecision(prec) << z1
         << "," << setprecision(prec) << z2
         << "," << setprecision(prec) << vx1
         << "," << setprecision(prec) << vx2
         << "," << setprecision(prec) << vy1
         << "," << setprecision(prec) << vy2
         << "," << setprecision(prec) << vz1
         << "," << setprecision(prec) << vz2
         << endl;
  }

  return 0;
}
