#ifndef __penningtrap_hpp__
#define __penningtrap_hpp__

#include <armadillo>
#include <vector>
#include <string>

// Including the particle header-file
#include "particle.hpp"

class PenningTrap{

// Make everything public so the code doesn't complain:
public:
	// Memeber variables
	double B0;	// Magnetic field strength
	double V0;	// Applied potential
	double d;	// Characteristic dimention
	std::vector<Particle>particles;	// Container of Particle objects

	// Constructor:
	PenningTrap(double B0_in, double V0_in, double d_in);

	// Defining member functions:
	// Adding a particle to the trap
	void add_particle(Particle p_in);

	// External electric field
	arma::mat external_E_field(arma::mat R);

	// External magnetic field
	arma::mat external_B_field(arma::mat V);

	// Force due to the interaction among the particles, particle_i and particle_j
	arma::vec force_particle(arma::vec r_i, arma::vec r_j, double q_i, double q_j);

	// Total force on particle_i from external fields
	arma::mat total_force_external(arma::mat R, arma::mat V);

	// Total force on particle from other particles
	arma::mat total_force_particles(arma::mat R);

	// Total force on particle_i from other particles and the external fields
	arma::mat total_force(arma::mat R, arma::mat V, bool particle_interaction);

	// Evolving the system with the Runge-Kutta 4th order
	void evolve_RK4(double dt, bool particle_interaction);

	// Evolving the system with forward Euler
	void forward_euler(double dt, bool particle_interaction);

};

#endif
