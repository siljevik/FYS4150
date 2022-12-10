//////////////////////////////////////////////////////////////////////////////////////
//                                                                                  //
// Compile with: g++ -O2 main.cpp -std=c++11 -I include -o main.exe -larmadillo     //
//                                                                                  //
// Run with: ./main.exe                                                             //
//                                                                                  //
//////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <armadillo>
#include <vector>

#include "funcs.hpp"
#include "funcs.cpp"

using namespace std;

int main(){

	funcs funcs;

	/*====================================*/
    /*~~~~ Constants, Variables, etc. ~~~~*/
    /*====================================*/
	int M = 5; // M-2 = 3 Here
	double h = 1;
	double dt = 1;
	double r_val = 1.0;
	int L = pow(M-2,2);

	/*====================================================*/
    /*~~~~ Making matrices filled with ones and zeros ~~~~*/
    /*====================================================*/
	arma::mat V(M, M, arma::fill::ones);
	arma::mat A(L, L, arma::fill::zeros);
	arma::mat B(L, L, arma::fill::zeros);

	// Define an empty vector that goes through a function that returns
	// a full vector dependent on indices (i,j) that we can run through with a loop
	cout << "Testing with matrix filled with ones for problem 2:" << endl;
	//cout << head.index_translator(M) << endl;



	////////////// SILJE TESTYTESTS //////////////
	//cout << "A:\n" << A;
	//////////////////////////////////////////////
	//head.matrix_filler(M,r_val);
	funcs.matrix_filler( M, r_val, L, A, B);
	funcs.diagonal_fill_AB(M, h, dt, L, V, A, B);
	arma::vec b = funcs.Bu_b(M, L, V, B);
 	return 0;
}
