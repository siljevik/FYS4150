#include <cmath>

// fetch class definition
#include "../header/particle.hpp"
#include "../header/penningtrap_t.hpp"

        // Constructor
        PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in, double f_in, double w_V_in){
	// define member variables as inputvariables
	B0 = B0_in;
	V0 = V0_in;
	d = d_in;
	f = f_in;
	w_V = w_V_in;
	}

        // Add a particle to the trap
        void PenningTrap::add_particle(Particle p_in){
		particles.push_back(p_in); // c++ append function to add particles
	}

        // External electric field at point r=(x,y,z) (E = -nabla V)
        arma::mat PenningTrap::external_E_field(double t, arma::mat R){
	int nn = particles.size(); //size of particle vector
	arma::mat E = arma::mat(3,nn);	// El.field strength matrix
	double  V0_t = V0*(1+f*cos(w_V*t)); // timedependent V0
	const double tr = V0 / (d*d); // timedependant rati ratio (tr) where V0 appear, should return 9.65 [u/(mu*u)**2 *e]
	// loop through positionvector R
	for (int i = 0; i < nn; i++){
	arma::vec r = R.col(i);
	if(arma::norm(r)){
	E.col(i) = -tr * arma::vec{-r(0) , -r(1), 2.* r(2)};
	}

	else{
	E.col(i) = arma::vec{0.0, 0.0, 0.0};
	}
	}
	return E;
	}

        // External magnetic field at point r=(x,y,z)
        arma::mat PenningTrap::external_B_field(arma::mat R, arma::mat V){
	int nn = particles.size();
	arma::vec B = arma::vec{0.0,0.0,B0}; // magnetic fieldstrength-vector
	arma::mat B_tot = arma::mat(3,nn); // total magnetic field-matrix


	// loop through the velocityvector V and update the charge and B_tot matrix
	for(int i=0; i<nn;i++){

	arma::vec r_i = R.col(i);

	if(arma::norm(r_i) < d){
	arma::vec v_i = V.col(i);
	double q_i = particles[i].q; // charge of particle i

	B_tot.col(i) = arma::cross(q_i*v_i, B); // qvxB
	}
	else{
	B_tot.col(i) = arma::vec{0.0, 0.0, 0.0};
	}
	}

	return B_tot;
	}

        // Force on particle_i from particle_j
        arma::vec PenningTrap::force_particle(arma::vec r_i, arma::vec r_j, double q_i, double q_j){
	double const ke = 1.38935333*pow(10,5); // Coulombs const.

	// calculate the electric field strength between the particles, with eq.(7)
	// from the project description

	double r_norm = arma::norm(r_i - r_j,1);
	arma::vec E = ke* q_j*(r_i-r_j)/pow(r_norm, 3);

	// returnig the force on particle i from particle j
	return (q_i * E);
	}


        // The total force on particle_i from the external fields (eq.(9) from project)
        arma::mat PenningTrap::total_force_external(double t, arma::mat R, arma::mat V){
	// Want to implement the Lorentz force on a particle from external fields
	// (F = qE + qvxB) by utilizing the external field functions for E and B:
	arma::mat extE = external_E_field(t, R); // E-field matrix
	arma::mat extB = external_B_field(R, V); // B-field matrix

	arma::mat ext_F_tot = extE + extB;// Total force matix

	return ext_F_tot;
	}

        // The total force on particle_i from the other particles
        arma::mat PenningTrap::total_force_particles(arma::mat R){

	int nn = particles.size(); //size of particle vector

	arma::mat F_tot = arma::mat(3,nn); // total force matrix
	// loop throug
	for(int i=0;i<nn;i++){
	arma::vec r_i = R.col(i);
	double q_i = particles[i].q;

		for(int j=0;j<nn;j++){
		if(j==i) {}
		else{
		double q_j = particles[j].q;
		arma::vec r_j = R.col(j);

		F_tot.col(i) += PenningTrap::force_particle(r_i, r_j, q_i, q_j);
		}
		}
	}
	return F_tot;
	}

        // The total force on particle_i from both external fields and other particles
        arma::mat PenningTrap::total_force(double t, arma::mat R, arma::mat V,bool particle_interaction){

	arma::mat F_tot; // total force matrix

	if(particle_interaction){
	F_tot = total_force_external(t, R, V) + total_force_particles(R);
	}
	else{
	F_tot = total_force_external(t, R,V);
	}
	return F_tot;
	}



	// Problem 7: Forward Euler and Runge-Kutta

        // Evolve the system one time step (dt) using Forward Euler
        void PenningTrap::evolve_forward_Euler(double t, double dt, bool particle_interaction){
	int nn = particles.size();
	arma::mat R = arma::zeros(3, nn);
	arma::mat V = arma::zeros(3, nn);
	arma::mat a = arma::zeros(3, nn);
	double m = particles[0].m;

	for(int i=0;i<nn;i++){
	R.col(i) = particles[i].r;
	V.col(i) = particles[i].v;
	}

	a = (1/m)*( total_force(t, R, V, particle_interaction));

	R += dt*V;
	V += dt*a;

	for(int i=0;i<nn;i++){
	particles[i].r = R.col(i);
	particles[i].v = V.col(i);
	}
}

        // Evolve the system one time step (dt) using Runge-Kutta 4th order
        void PenningTrap::evolve_RK4(double t, double dt, bool particle_interaction){
	int nn = particles.size();

	arma::mat R = arma::zeros(3, nn); // directional matrix
        arma::mat V = arma::zeros(3, nn); // velocity matrix
        arma::mat a = arma::zeros(3, nn); // a matrix
	double m = particles[0].m;

	arma::mat K1_r, K2_r, K3_r, K4_r; // Kmatrices for direction
	arma::mat K1_v, K2_v, K3_v, K4_v; // Kmatrices for velocity

	for(int i=0;i<nn;i++){
	R.col(i) = particles[i].r;
	V.col(i) = particles[i].v;
	}
	// K1
	a = (1/m)*( total_force(t,R, V, particle_interaction) );
	K1_v = dt*a;
	K1_r = dt*V;
	// K2
	a = (1/m)*( total_force(t + 0.5*dt, R+0.5*K1_r, V+0.5*K1_v, particle_interaction) );
        K2_v = dt*a;
        K2_r = dt*(V+0.5*K1_v);
	// K3
	a = (1/m)*( total_force(t + 0.5*dt, R+0.5*K2_r, V+0.5*K2_v, particle_interaction) );
        K3_v = dt*a;
        K3_r = dt*(V+0.5*K2_v);
	// K4
	a = (1/m)*( total_force(t + dt, R+0.5*K3_r, V+0.5*K3_v, particle_interaction) );
        K4_v = dt*a;
        K4_r = dt*(V+0.5*K3_v);
	// updating the R and V matrices
	R += (1.0/6)*(K1_r + 2*K2_r + 2*K3_r + K4_r);
	V += (1.0/6)*(K1_v + 2*K2_v + 2*K3_v + K4_v);

	for(int i=0;i<nn;i++){
	particles[i].r = R.col(i);
	particles[i].v = V.col(i);
	}
	}


// adding a function that counts the particles in the trap
	int PenningTrap::count_particles()
	{
	int n = particles.size();
	return n;
	}

        // all of the memnber functions above will be implementet in the functionfile
	// ../func/penningtrap.cpp

