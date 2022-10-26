//////////////////////////////////
//
// Compile with: g++ -O2 main.cpp func/particle.cpp -std=c++11 -I include -o main.exe -larmadillo
//
// Run with: ./main.exe
//
/////////////////////////////////
#include <iostream>
#include <armadillo>
#include <cmath>
#include <chrono>
#include <fstream>
#include <iomanip>


// namespaces
using namespace std;

// including the header files containing our classes
#include "header/particle.hpp"
#include "header/penningtrap.hpp"
// including our file containing all the functions
#include "func/penningtrap_t.cpp"

using namespace std;


int main(){
	// constants:
	double const ke = 1.38935333*pow(10,5); // Coloumb constant [(u*(mu*m)e3) / (mu*s)e2 * ee2 ]
	// derived SI units of:
	double const T = 9.64852558*pow(10, 1); // magnetic field strength Tesla
	double const V = 9.64852558*pow(10, 7); // electric potential Volt
	// default Penning Trap configurations:
	double const B0 = 9.65*pow(10,1); // default mag. field strength of penningtrap [u/((mu*s)*e)]
	double const V0 = 2.41*pow(10, 6); // electric potential [((u*(mu*m)e2)/(mu*s)e2 *e)]
	double const d = 500;		// mu*m
	double const mCa_pos = 40.078;	// mass of Ca+ ion [u]
	double const q = 1;		// charge of +1 e (elementary charge)

/////////////// problem 9

	arma::vec f = arma::vec{0.1, 0.4, 0.7}; // amplitudes
	double w_V = 0.2;	// angular frequency [MHz]

	PenningTrap trap = PenningTrap(B0, V0, d, f(0), w_V);

	// Note; to set the seed for Armadillos random number generator:
	arma::arma_rng::set_seed_random();

	int t_tot = 500; // microseconds
	double dt = 1.3*pow(10, -2);
	double t_steps = ceil(t_tot/dt);

	double f_max = 2.5;
	double f_min = 0.2;
	double dw = 0.02;
	double f_steps = ceil((f_max - f_min)/dw);

	arma::vec t = arma::linspace(0, t_tot, t_steps);
	arma::vec f_list = arma::linspace(f_min, f_max, f_steps);
	arma::mat n_particles = arma::mat(f_steps, 3);



	// running simulation
//all amplitudes empties the trap, using only f=0.1
for(int i=0; i<1; i++)
{
	//setting amplitude
	trap.f = f(i);

	// timing algorithm
	auto t1 = chrono::high_resolution_clock::now();
		for(int j=0; j < f_steps; j++)
		{
		//removing all particles from the trap
		trap.particles.clear();

			// filling the trap
			for(int k=0; k<50; k++)
			{
			arma::vec r = arma::vec(3).randn() * 0.1 * trap.d;
			arma::vec v = arma::vec(3).randn() * 0.1 * trap.d;

			Particle p = Particle(q, mCa_pos, r, v); //paricle(Ca+) properties
			trap.add_particle(p);
			}


		//staring with 100 particles
		n_particles.col(i)(0) = trap.particles.size();

		// setting the frequency
		trap.w_V = f_list(j);


			for(int l = 0; l< t_steps; l++)
			{
			trap.evolve_RK4(t(l), dt, false);
			}
			if(j != 0)
			{
			n_particles.col(i)(j) = trap.count_particles();
			}
		}

		auto t2 = chrono::high_resolution_clock::now();
		double duration = chrono::duration<double>(t2-t1).count();

		// printing progress to screen
		cout 	<< "Finished " << i+1 << "/3" << "in"
			<< int(duration) / 60 << " min, "
			<< int(duration) % 60 << " s" << endl;
}

/*
	// getting the fraction of remeaining particles
	n_particles /= 50;

	int width = 15;
	int prec = 6;

	ofstream datafile;
	datafile.open("no_of_particles.txt", ofstream::out | ofstream::trunc);

	datafile << "#" << setw(width - 1) << "freq"
                 << setw(width) << "f1"
                 << setw(width) << "f2"
                 << setw(width) << "f3"
                 << setw(width) << "n1 left"
                 << setw(width) << "n2 left"
                 << setw(width) << "n3 left"
                 << endl;
	for (int i = 0; i < f_steps; i++)
	{
	    datafile << setw(width) << setprecision(prec) << f_list(i)
            << setw(width) << setprecision(prec) << f(0)
            << setw(width) << setprecision(prec) << f(1)
            << setw(width) << setprecision(prec) << f(2)
            << setw(width) << setprecision(prec) << n_particles.col(0)(i)
            << setw(width) << setprecision(prec) << n_particles.col(1)(i)
            << setw(width) << setprecision(prec) << n_particles.col(2)(i)
            << endl;
  	}
 	 datafile.close();
*/
return 0;
}
