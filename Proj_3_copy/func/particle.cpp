// fetch the particle class from the hpp-file
#include "../header/particle.hpp"

// constructor
Particle::Particle(double chrg_in, double mass_in, arma::vec pos_in, arma::vec v_in){
	// assigning the member variables
	q = chrg_in;
	m = mass_in;
	r = pos_in;
	v = v_in;
}

// To write out the information about the particle in a string
std::string Particle::p_info(){
std::string p_props = "Properties of particle: \nCharge " + std::to_string(q) + " e.\nMass "
			+ std::to_string(m) + " u";

return p_props;
}



