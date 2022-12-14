//////////////////////////////////
//
// Compile with: g++ single_test.cpp func/particle.cpp func/penningtrap.cpp func/analytical.cpp -std=c++11 -I include -o single_test.exe -larmadillo
//
// To execute file and save data: ./single_test.exe > single_particle_z.txt
//
/////////////////////////////////
#include <iostream>
#include <cmath>
#include <armadillo>
#include <iomanip> // needed for the setw() and setprecision()


// class definition
#include "header/particle.hpp"
#include "header/penningtrap.hpp"
#include "header/analytical.hpp"

//////////////////////////////
// ~~~~~~~~~~~~~~~~~~~~~~~~ //
// TEST FOR SINGLE PARTICLE //
// ~~~~~~~~~~~~~~~~~~~~~~~~ //
//////////////////////////////
using namespace std;

int main()
{
  // int width = 18; // We don't use it ?
  int prec = 4;

  PenningTrap trap = PenningTrap(B0, V0, d);

  double q = 1.0;     // Charge of particle [e]
  double m = 40.078;  // Mass of Ca+ ion [u]

  // Initial values of Particle 1
  double x_0 = 20.0, y_0 = 0.0, z_0 = 20.0;
  double vx_0 = 0.0, vy_0 = 25.0, vz_0 = 0.0;

  arma::vec r = arma::vec{x_0, y_0, z_0};
  arma::vec v = arma::vec{vx_0, vy_0, vz_0};

  Particle p = Particle(q, m, r, v);

  trap.add_particle(p);

<<<<<<< HEAD
  /*int t_tot = 50;         // Total time
  double dt = pow(10,-3);      // Stepsize time
  int steps = t_tot / dt; // Number of steps*/

  int t_tot = 50;         // Total time
  int steps = 1000; // Number of steps
  double dt = t_tot/steps;      // Stepsize time
  
=======
<<<<<<< HEAD
  int t_tot = 50; // microseconds
=======
  int t_tot = 5000;
>>>>>>> cf16f5cc0465c184ee1715ed0f11d9846509ad83
  double dt = pow(10,-3);
  int steps = t_tot / dt;
>>>>>>> 5b78f0612ba4a400a575a89553b4856d1f507187

  arma::vec t = arma::linspace(0, t_tot, steps);

  // num = numerical
  // ana = analytical
  cout << "Time"
       << "," << "x_num"
       << "," << "y_num"
       << "," << "z_num"
       << "," << "x_ana"
       << "," << "y_ana"
       << "," << "z_ana"
       << endl;
     // Printing the values with a comma (",") between them
  cout << setprecision(prec) << t(0)
       << "," << setprecision(prec) << x_0
       << "," << setprecision(prec) << y_0
       << "," << setprecision(prec) << z_0
       << "," << setprecision(prec) << x_0
       << "," << setprecision(prec) << y_0
       << "," << setprecision(prec) << z_0
       << endl;

	for (int i = 1; i < steps-1; i++)
  	{
    	trap.evolve_RK4(dt, true);

    	arma::vec r_ana = analytical_solution(r, v, p.q, B0, p.m, V0, d, t(i));

    	double x_ana = r_ana[0], y_ana = r_ana[1], z_ana = r_ana[2];
    	double x_num = trap.particles[0].r[0];
    	double y_num = trap.particles[0].r[1];
    	double z_num = trap.particles[0].r[2];

    	cout << setprecision(prec) << t(i)
         << "," << setprecision(prec) << x_num
         << "," << setprecision(prec) << y_num
         << "," << setprecision(prec) << z_num
         << "," << setprecision(prec) << x_ana
         << "," << setprecision(prec) << y_ana
         << "," << setprecision(prec) << z_ana
         << endl;
  	}
  return 0;
}
