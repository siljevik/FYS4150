#ifndef __penningtrap_hpp__
#define __penningtrap_hpp__

#include <armadillo>
#include <string>
#include <vector>
// include the particle class-file
#include "particle.hpp"


//this is the pennigtrap class
class PenningTrap{
	public:

	double B0;	// Magnetic field strangth
	double V0;	// applied potential
	double d;	// characteristic dimesnion
	double f;	// fr
	double w_V;
	std::vector<Particle> particles;	//vector that contains all the particle objects


	// Adding the suggested codesnippet for the PenningTrap memberfunctions, from the project text (for the ease of mind):

	// Constructor
	PenningTrap(double B0_in, double V0_in, double d_in, double f_in, double w_V_in);

	// Add a particle to the trap
	void add_particle(Particle p_in);

	//////////////////// Problem 6 member functions
	// External electric field at point r=(x,y,z)
	arma::mat external_E_field(double t, arma::mat R);

	// External magnetic field at point r=(x,y,z)
	arma::mat external_B_field(arma::mat R, arma::mat V);

	// Force on particle_i from particle_j
	arma::vec force_particle(arma::vec r_i, arma::vec r_j, double q_i, double q_j);

	//////////////////// Problem 7 member function
	// The total force on particle_i from the external fields
	arma::mat total_force_external(double t, arma::mat R, arma::mat V);

	// The total force on particle_i from the other particles
	arma::mat total_force_particles(arma::mat R);

	// The total force on particle_i from both external fields and other particles
	arma::mat total_force(double t, arma::mat R, arma::mat V, bool particle_interaction);

	// Evolve the system one time step (dt) using Runge-Kutta 4th order
	void evolve_RK4(double t,double dt, bool particle_interaction);

	// Evolve the system one time step (dt) using Forward Euler
	void evolve_forward_Euler(double t, double dt, bool particle_interaction);


	// adding a member function that counts the number of particles
	int count_particles();
	// all of the memnber functions above will be implementet in the functionfile ../func/penningtrap.cpp
};
#endif
