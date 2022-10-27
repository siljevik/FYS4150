//////////////////////////////////
//
// Compile with: g++ r_err.cpp func/particle.cpp func/penningtrap.cpp func/analytical.cpp -std=c++11 -I include -o r_err.exe -larmadillo
//
/////////////////////////////////
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip> // needed for the setw() and setprecision()

// classfiles
#include "header/particle.hpp"
#include "header/penningtrap.hpp"
#include "header/analytical.hpp"

using namespace std;

// main function of the r_err
int main()
{
  int width = 15;
  int prec = 4;

  PenningTrap trap = PenningTrap(B0, V0, d);

  double q = 1.0; // charge  
  double m = 40.078;  // mass of Ca+ ion in u

  double x_0 = 1.0, y_0 = 0.0, z_0 = 1.0;
  double vx_0 = 0.0, vy_0 = 1.0, vz_0 = 0.0;

  arma::vec r = arma::vec{x_0, y_0, z_0};
  arma::vec v = arma::vec{vx_0, vy_0, vz_0};

  Particle p = Particle(q, m, r, v);

  double t_tot = 50;
  arma::vec dMax = arma::vec(4);
  arma::vec h = arma::vec(4);

  ofstream outfile;
  outfile.open("data_r_err_FE.txt", ofstream::out | ofstream::trunc);

  outfile << "#" << "," << "hk"
          << "," << "time"
          << "," << "r_ana"
          << "," << "r_num"
          << "," << "r_err"
          << endl;

  for (int i = 1; i < 5; i++)
  {
        // run for nk = 4000, 8000, 16000, 32000
    int nk = 2000*pow(2,i);        //int steps = 2000 * pow(2, i);
    double hk = t_tot/nk; //double dt = t_tot / steps;

    h(i-1) = hk;

    arma::vec t = arma::linspace(0, t_tot, nk);

    trap.particles.clear();
    trap.add_particle(p);
    trap.particles[0].r = r;

    outfile << "," << setprecision(prec) << hk
            << "," << setprecision(prec) << t(0)
            << "," << setprecision(prec) << arma::norm(r)
            << "," << setprecision(prec) << arma::norm(r)
            << "," << setprecision(prec) << 0
            << endl;

    for (int j = 1; j < nk; j++)
    {

      // run this for RK4 and forward_euler
      trap.forward_euler(hk, true);
      // trap.evolve_RK4(hk, true);

      // Numerical and analytical position of particle 1
      arma::vec r_num = trap.particles[0].r;
      arma::vec r_ana = analytical_solution(r, v, p.q, B0, p.m, V0, d, t(j));

      // Absolut and relative error
      double abs_err = arma::norm(r_num - r_ana);
      double rel_err = abs_err / arma::norm(r_ana);
      
      if (abs_err > dMax(i - 1))
      {
        dMax(i - 1) = abs_err;
      }

      outfile << "," << setprecision(prec) << hk
              << "," << setprecision(prec) << t(j)
              << "," << setprecision(prec) << arma::norm(r_ana)
              << "," << setprecision(prec) << arma::norm(r_num)
              << "," << setprecision(prec) << rel_err
              << endl;
    }
  }
  outfile.close();


////////////////////////////
//                        //
//  Last point of task 8  //
//                        //
////////////////////////////

// Want to use the data from the simulation to estimate the convergence rate of RK4- and FE-methods
  double r_err = 0.;


// Run through the simulation results and calculate with the relative error equation from the project description
  for (int k = 1; k < 4; k++)
  {
    r_err += (1.0/3) * (log(dMax(k) / dMax(k-1)) / log(h(k) / h(k-1)));
  }

// Print out the convergence rate in the terminal:
   cout << "Convergence rate for forward Euler: " << setprecision(4) << r_err << endl;
 // cout << "Convergence rate for RK4: " << setprecision(4) << r_err << endl;

  /* From terminal:
  > Convergence rate for FE: 1.443
  > Convergence rate for RK4: 1
  */

  return 0;
}
