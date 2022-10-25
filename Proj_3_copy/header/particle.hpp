#ifndef __particle_hpp__
#define __particle_hpp__

#include <armadillo>
#include <string>


// class of particle where we define member variables that will be our particle properties
class Particle {
public:
	double q;	//charge of particle [e]
	double m;	//mass [u]
	arma::vec r;	// position-vector [micrometre]
	arma::vec v;	// velocity-vector [m/s]

// The constructor should assign values to the member variables
	Particle(double q_in, double m_in, arma::vec r_in, arma::vec v_in); //constructor of member variables

	std::string p_info();
};

#endif
