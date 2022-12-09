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

	private: // Nothing to see here.

	public:

	arma::vec vectorfiller(int M);

	void index_translator(int M, int k, int & i, int & j);


	void matrixfiller(int M, std::complex<double> r_val);

}; // end of class Header

#endif
