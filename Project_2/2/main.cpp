#include <iostream>
#include <armadillo>
// wirite a code that set up the tridiagonal, 6X6 matrix A.
// solves Av = hv using armadillo's arma::eig_sym,
// checks that the eigenvalues and eigenvectors from Armadillo agrees with the analytical
// results for N = 6

using namespace arma;

int main() {
// big N to adjust the size of the NxN matrix A
	int N = 6;
// minimum and maximum points
	double x_min = 0.0;
	double x_max = 1.0;

	int  n = N + 1;
	double dx = x_max - x_min;
	double h = dx*1./n;
// create an empty matrix A, using Armadillo
	mat A(N, N, fill::zeros);

// Creating two loops to fill the tridiagonal matrix with values:
// values, a ,diag.
for (int i=0; i<N; i++){
 A(i, i) =  2./(h*h);
}

// values,c ,diag.
for (int i=0; i<(N-1) ; i++){
 A(i+1, i) = -1./(h*h);
 A(i, i+1) = -1./(h*h);
}

// we want to solve Av = lambda v using arma::eig_sym

	vec eigval; // vector containing eigenvalues
	mat eigvec; // eigenvector
	eig_sym(eigval, eigvec, A); // gives the correct answer

// printing out the tridiagonal matrix
A.print(std::cout);
// eigenvalues in order
eigval.print(std::cout);
// eigevectors in a matrix
eigvec.print(std::cout);

return 0;
}
