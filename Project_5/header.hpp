#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <random>
#include <map>

#ifndef __header_hpp__
#define __header_hpp__

class Header
{

	private: 
	// Nothing to see here.

	public:

	arma::vec vector_filler(int M, arma::mat V);

	void index_translator(int M, int k, int & i, int & j);

	void matrix_filler(int M, double r_val);

	void diagonal_fill_AB(int M, double h, double dt, arma::mat V,arma::mat & A, arma::mat & B);

}; // end of class Header

#endif
