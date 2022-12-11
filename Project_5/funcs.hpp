#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <random>
#include <map>

#ifndef __funcs_hpp__
#define __funcs_hpp__

class funcs
{

	private: 
	// Nothing to see here.

	public:

	arma::vec vector_filler(int M, arma::mat V);

	void index_translator(int M, int k, int & i, int & j);

	void matrix_filler(int M, double r_val, int L, arma::mat & A, arma::mat & B);

	void diagonal_fill_AB(int M, double h, double dt, int L, arma::mat V,arma::mat & A, arma::mat & B);

	arma::vec Bu_b(int M, int L, arma::mat V, arma::mat B);

	arma::vec Au_b(arma::mat A, arma::vec b);

	void initial_u(int M, double h, int L, arma::vec u_0);

}; // end of class Header

#endif
