#include <iostream>
#include <armadillo>
#include <vector>

#include "header.hpp"
#include "funcs.cpp"

using namespace std;

int main(){

	Header head;

	int M = 5; // M-2
	std::complex<double> r_val = 1.0;
	// Define an empty vector that goes through a function that returns
	// a full vector dependent on indices (i,j) that we can run through with a loop
	cout << "Testing with matrix filled with ones for problem 2:" << endl;
	cout << head.indextranslator(M) << endl;

	head.matrixfiller(M,r_val);
 return 0;

}
