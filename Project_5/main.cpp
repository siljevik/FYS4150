#include <iostream>
#include <armadillo>
#include <vector>

#include "header.hpp"
#include "funcs.cpp"

using namespace std;

int main(){

	Header head;

	/*====================================*/
    /*~~~~ Constants, Variables, etc. ~~~~*/
    /*====================================*/
	int M = 5; // M-2 = 3 Here
	int h = 1;
	int dt = 1;
	arma::cx_double r_val = 1.0;

	/*===========================================*/
    /*~~~~ Making matrices filled with zeros ~~~~*/
    /*===========================================*/
	arma::mat V(M,M,fill::ones);
	arma::mat A(pow(M-2,2), pow(M-2,2), fill::zeros);
	arma::mat B(pow(M-2,2), pow(M-2,2), fill::zeros);

	// Define an empty vector that goes through a function that returns
	// a full vector dependent on indices (i,j) that we can run through with a loop
//	cout << "Testing with matrix filled with ones for problem 2:" << endl;
//	cout << head.index_translator(M) << endl;

	head.matrix_filler(M,r_val);

	diagonal_fill_AB(M, h, dt, V, & A, & B);

 	return 0;

}
