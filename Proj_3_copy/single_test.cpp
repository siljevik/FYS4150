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


// Test for a single particle

using namespace std;

int main()
{
  int width = 18;
  int prec = 4;

  PenningTrap trap = PenningTrap(B0, V0, d);

  double q = 1.0; // charge of particle [e]
  double m = 40.078;  //Mass of Ca+ ion [u]

  double x_0 = 20., y_0 = 0., z_0 = 20.;
  double vx_0 = 0., vy_0 = 25., vz_0 = 0.;

  arma::vec r = arma::vec{x_0, y_0, z_0};
  arma::vec v = arma::vec{vx_0, vy_0, vz_0};

  Particle p = Particle(q, m, r, v);

  trap.add_particle(p);

  int t_tot = 50; // microseconds
  double dt = pow(10,-3);
  int steps = t_tot / dt;

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
