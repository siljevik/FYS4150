#ifndef __penningtrap_hpp__
#define __penningtrap_hpp__

#include <armadillo>
#include <vector>
#include <string>

//include header
#include "particle.hpp"

class PenningTrap{

// make everything public so the code doesn't complain:
public:
	// memeber variables
	double B0;	// Magnetic field strength
	double V0;	// applied potential
	double d;	// characteristic dimention
	std::vector<Particle>particles;	// container of Particle objects

//constructor:
	PenningTrap(double B0_in, double V0_in, double d_in);

//defining member functions:
	// adding a particle to the trap
	void add_particle(Particle p_in);

	// external electric field
	arma::mat external_E_field(arma::mat R);

	// external magnetic field
	arma::mat external_B_field(arma::mat V);

	// force due to the interaction among the particles, particle_i and particle_j
	arma::vec force_particle(arma::vec r_i, arma::vec r_j, double q_i, double q_j);

	// total force on particle_i from external fields
	arma::mat total_force_external(arma::mat R, arma::mat V);

	// total force on particle from other particles
	arma::mat total_force_particles(arma::mat R);

	// total force on particle_i from other particles and the external fields
	arma::mat total_force(arma::mat R, arma::mat V, bool particle_interaction);

	// evolving the system with the Runge-Kutta 4th order
	void evolve_RK4(double dt, bool particle_interaction);

	// evolving the system with forward euler
	void forward_euler(double dt, bool particle_interaction);

};

#endif
