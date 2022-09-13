#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <ctgmath>

// main function
int main() {


	// conditions
	double x_min = 0.0;
	double x_max = 1.0;
	// steps/stepsize
	double n_steps = 100;
	double h = (x_max - x_min)/n_steps;

// we want to implement our specialized algorithm in code
// we will need the parametres v, g_tilde ,b_tilde, b and g
// in loops to run through the values needed for calculation

	// forward-substitution-loop
	for (int i=0; i<= n_steps; i++)
	{
	// x and f(x) used in calculation
	double x = x_min;
	x += h;
	double f = 100*exp(-10*x);

	//defining  g_tilde and b_tilde for later use
	double hh = h*h;
	double b_tilde(0) = 2;
	double b_tilde(i) = 2- 1/b_tilde(i-1);

	double g_tilde(0) = hh*f( x(0) );
	double g_tilde(i) = hh*(f( x(i) ) + 0.5*f( x(i-1) ));
	}

	// backwards-substitution

	double v(n_steps) = (g_tilde(n_steps) + v(n_steps-1))/(2-(1/b_tilde(n_steps-1)));
	for (int i=n_steps - 1; i>= 0; i--)
	{	//defining our specialized algorithm for counting backwards
		double v(i) = (g_tilde(i) + v(i-1))/(2- (1/b_tilde(i-1));
	}


	return 0;
}
