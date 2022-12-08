#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <random>
#include <vector>
#include <map>


#include "header.hpp"


arma::vec Header::indextranslator(int M){

        // defining an empty vector
        std::vector<double> un_vec;
	// test values
	int L = pow( (M-2), 2);
	std::cout << "Length: " << L << std::endl;
	arma::mat uij(L,L,arma::fill::ones);
	std::cout << uij << std::endl;
        // fill the vector with vectors
                for(int i = 0; i < M-2; i++)
                {       int column = i;
			for(int j = 0; j < M-2; j++)
                        {
                                //append uij values into the uij_vec vector which will be appended into un_vec
                                un_vec.push_back(uij(column,j));
                        }//end of j-loop
                }// end of i-loop
        return un_vec;
}

void Header::matrixfiller(int M, std::complex<double> r_val){

	// For making the matrices we want to use sp_cx_mat, as it is ideal for storing huge matrices containing mostly zeroes.
	int L = pow( (M-2), 2);
	// start out with two empty matrices A and B
	arma::sp_cx_mat A(L,L);
	arma::sp_cx_mat B(L,L);

	// Defining vector containing value r
	arma::cx_vec r_diagonal( (M-2)*(M-3) );// M-3 because we remove 1 position to shift the position in A/B.
	r_diagonal.fill( r_val ); // fill r_vec with r

	// Now we want to create the diagonals containing r (or -r) values
	arma::sp_cx_mat subdiagonal(L,L);
	arma::sp_cx_mat superdiagonal(L,L);

	subdiagonal.diag(-1);
	superdiagonal.diag(1);


	// want to insert the diaginal matrices containing r into A and B
	A.diag(M-2) = -r_diagonal;
	A.diag(2-M) = -r_diagonal;
        B.diag(M-2) = +r_diagonal;
        B.diag(2-M) = +r_diagonal;



	// Defining two vectos with cx_vec for the sake of complex numbers.
	arma::cx_vec a_vec(L, arma::fill::zeros);
	arma::cx_vec b_vec(L, arma::fill::zeros);

	a_vec.fill(1.0);
	b_vec.fill(1.0);
	// Adding the a and b vectors to the A and B matrices
//	A.diag(M-2) = a_vec;
//	B.diag(M-2) = b_vec;


	std::cout << "Matrix A " << "\n";
	std::cout << A << std::endl;
//	std::cout << B << std::endl;
}
