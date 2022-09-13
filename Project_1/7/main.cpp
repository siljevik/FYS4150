#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <ctgmath>

int main() {

// parametres(x-boundaries)
	double x_min = 0.0;
	double x_max = 1.0;
// step size
	int n_steps = 100;
	double h = (x_max - x_min)/n_steps;
// "arrays"
	// since the values on the diagonal is constant, we write the
	// arrays as constants to use in calculations of b_tilde, g_tilde and v
	int a = -1;
	int b =  2;
	int c = -1;
// functions
	//double x = 0;
	double hh = h*h;
	//double f = 100*exp(-10.*x);
// vectors
	double x(0) =  h;
	double b_tilde(0) =  b;
	double g_tilde(0) = hh*f( x[0] );

// forward substitution
	for (int i=1; i<= n_steps ; i++) {
	//filling arrays with values
	double x(i) = x(i-1) +h;
	double g(i) = hh*100*exp( -10*x(0) );

	// forward subst. part
	double b_tilde(i) = b - (a/b_tilde(i-1))*c;
	double g_tilde(i) = g(i) -(a/b_tilde(i-1))*g_tilde(i-1);
	}

// backwards substitution
	double v(n_steps-1) = g_tilde(n_steps - 1)/b_tilde(n_steps -1);
	for (int i = n_steps -2; i>=0; i--){
	double v(i) = (g_tilde(i) - v(i+1)*c)/b_tilde(i);
	}
	return 0;
}
