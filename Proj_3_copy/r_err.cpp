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

  double x0 = 1.0, y0 = 0.0, z0 = 1.0;
  double vx_0 = 0.0, vy_0 = 1.0, vz_0 = 0.0;

  arma::vec r = arma::vec{x0, y0, z0};
  arma::vec v = arma::vec{vx_0, vy_0, vz_0};

  Particle p = Particle(q, m, r, v);

  double t_tot = 50; // microseconds
  arma::vec dMax = arma::vec(4);
  arma::vec h = arma::vec(4);

  ofstream outfile;
  outfile.open("data_r_err.txt", ofstream::out | ofstream::trunc);

  outfile << "#" << setw(width) << "dt"
          << setw(width) << "time"
          << setw(width) << "r_ana"
          << setw(width) << "r_num"
          << setw(width) << "r_err"
          << endl;

  for (int i = 1; i < 5; i++)
  {
    // nk can be written as n_k = 2000 * 2^{k}, k = 1,2,3,4
    int steps = 2000 * pow(2, i);
    double dt = t_tot / steps;

    h(i-1) = dt;

    arma::vec t = arma::linspace(0, t_tot, steps);

    trap.particles.clear();
    trap.add_particle(p);
    trap.particles[0].r = r;

    outfile << setw(width) << setprecision(prec) << dt
            << setw(width) << setprecision(prec) << t(0)
            << setw(width) << setprecision(prec) << arma::norm(r)
            << setw(width) << setprecision(prec) << arma::norm(r)
            << setw(width) << setprecision(prec) << 0
            << endl;

    for (int j = 1; j < steps; j++)
    {

      // run this for RK4 and forward_euler
    //   trap.forward_euler(dt, true);
    trap.evolve_RK4(dt, true);

      arma::vec r_num = trap.particles[0].r;
      arma::vec r_ana = analytical_solution(r, v, p.q, B0, p.m, V0, d, t(j));

      double abs_err = arma::norm(r_num - r_ana);
      double rel_err = abs_err / arma::norm(r_ana);

      if (abs_err > dMax(i - 1))
      {
        dMax(i - 1) = abs_err;
      }

      outfile << setw(width) << setprecision(prec) << dt
              << setw(width) << setprecision(prec) << t(j)
              << setw(width) << setprecision(prec) << arma::norm(r_ana)
              << setw(width) << setprecision(prec) << arma::norm(r_num)
              << setw(width) << setprecision(prec) << rel_err
              << endl;
    }
  }
  outfile.close();

  double r_err = 0.;

  for (int k = 1; k < 4; k++)
  {
    r_err += 1. / 3 * (log(dMax(k) / dMax(k-1)) / log(h(k) / h(k-1)));
  }

  // cout << "Convergence rate for forward Euler: " << setprecision(4) << r_err << endl;
   cout << "Convergence rate for RK4: " << setprecision(4) << r_err << endl;

  /* From terminal:
  > Convergence rate for FE: 1.395
  > Convergence rate for RK4: 1
  */

  return 0;
}
