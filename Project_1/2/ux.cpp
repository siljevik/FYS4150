#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <ctgmath>

int main(){
// name of data-file
	std::string filename = "x_u.txt";
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
	int width = 12;
	int prec = 4;
//initial x and u values
	double x = x_min;
	double u = 1-(1-exp(-10))*x - exp(-10*x);
	//loop over steps
	for (int i=0; i<= n_steps; i++)
	 {
	//write a line with the current x and y values to file
	ofile << std::setw(width) << std::setprecision(prec) << std::scientific << x
	      << std::setw(width) << std::setprecision(prec) << std::scientific << u
	      << std:: endl;
	// update x and u(x) values
	x += h;
	u = 1-(1-exp(-10))*x - exp(-10*x);
	}

	//close the output file
	ofile.close();
	//exit program
	return 0;
}
