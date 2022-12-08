#include <iostream>
#include <armadillo>
#include <vector>

#include "header.hpp"
#include "funcs.cpp"

using namespace std;

int main(){

	Header head;

	int M = 3; // M-2
	// cx_mat
	arma::mat A(M,M,arma::fill::randn);
	arma::mat B = A + A;

	arma::cx_mat X(A,B);
	/*cout << A << endl;
	cout << B << endl;
	cout << X << endl;*/

	//cx_double
	X(1,2) = arma::cx_double(2.0, 3.0);

	cout << X << endl;
	cout << X(1,2) << endl;

	// Define an empty vector that goes through a function that returns
	// a full vector dependent on indices (i,j) that we can run through with a loop
	cout << head.indextranslator(M) << endl;
 return 0;

}
