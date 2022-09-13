#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <ctgmath>

int main(){
// name of data-file
	std::string filename = "x_v.txt";
	 // Create and open the output file. Or, technically, create
  // an "output file stream" (type std::ofstream) and connect it to our filename.
	std::ofstream ofile;
	ofile.open(filename);

	// Set some parameters for our computation
	double x_min = 0.0;
	double x_max = 1.0;
// step size
	int n_steps = 100;
	double h = (x_max - x_min) / n_steps;

// parameters to adjust the data and significant numbers in the txt file
	int width_x = 8;
	int width_v = 11;
	int prec = 4;
//initial x and u values
	double x = x_min;
	double u = (g-v*c)/b;
	//loop over steps
	for (int i=0; i<= n_steps; i++)
	 {
	//write a line with the current x and y values to file
	ofile << std::setw(width_x) << std::setprecision(prec) << std::scientific << x
	      << std::setw(width_v) << std::setprecision(prec) << std::scientific << v
	      << std:: endl;
	// update x and u(x) values
	x += h;
	u = (g-v*c)/b;
// It doesn't stop at x=1. WHY?!?!
	}

	//close the output file
	ofile.close();
	//exit program
	return 0;
}