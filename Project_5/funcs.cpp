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

using namespace std;

/*===================================*/
/*~~~~~~     Vector Filler     ~~~~~~*/
/*===================================*/
arma::vec Header::vector_filler(int M){
    // Defining an empty vector
    std::vector<double> un_vec;
	//////////// For code testing ////////////
	int L = pow( (M-2), 2);
	std::cout << "Length: " << L << std::endl;
	arma::mat uij(L,L,arma::fill::ones);
	std::cout << uij << std::endl;
	/////////////////////////////////////////

	// fill the vector with vectors
    for(int i = 0; i < M-2; i++){
		int column = i;
		for(int j = 0; j < M-2; j++){
			// Append uij values into the uij_vec vector which will be appended
			// into un_vec
            un_vec.push_back(uij(column,j));
        }//end of j-loop
    }// end of i-loop
    return un_vec; // Returns the vector
}



/*======================================*/
/*~~~~~~     Index Translator     ~~~~~~*/
/*======================================*/
void Header::index_translator(int M, int k, int & i, int & j){
	// Test
	int column_from_k = (k+1)/(M-2); // int/int will make 5/2 = 2, or 6/7 = 0
	i = k-(column_from_k*(M-2))+1;
	j = column_from_k + 1;
}



/*===================================*/
/*~~~~~~     Matrix Filler     ~~~~~~*/
/*===================================*/
void Header::matrix_filler(int M, double r_val){

	// For making the matrices we want to use sp_cx_mat, as it is ideal for storing huge matrices containing mostly zeroes.
	int L = pow( (M-2), 2);
	// start out with two empty matrices A and B
	arma::mat A(L,L, arma::fill::zeros);

	// Defining vector containing value r
	arma::vec r_outer_diagonal( (M-2)*(M-3), arma::fill::none );// M-3 because we remove 1 position to shift the position in A/B.
	r_outer_diagonal.fill( r_val );

	// Checkpoint
	std::cout << "Outer diagonal vector: " << "\n";
	std::cout << r_outer_diagonal << std::endl;

	// want to insert the diaginal matrices containing r into A and B
	A.diag(M-2) = -r_outer_diagonal;
	A.diag(2-M) = -r_outer_diagonal;

	// Checkpoint
	std::cout << "Adding outer diagonal to A: " << "\n";
	std::cout << A << std::endl;


	// inner diagonal vectors containng r and for each 3rd position,0.
        arma::vec r_inner_diagonal( (M-1)*(M-3), arma::fill::none );// M-3 because we remove 1 position to shift the position in A/B.
	r_inner_diagonal.fill(r_val);

        for(int k=0; k < (M-1)*(M-3); k++)
        {
                if( (k+1)%3 == 0 )
                {
                        r_inner_diagonal[k] = 0;
                } // end of if-statement
        }// end of k-loop

	// Checkpoint
	std::cout << "r_inner_diagonal vector: " << "\n" ;
	std::cout << r_inner_diagonal	<< std::endl;

	A.diag(1)  = -r_inner_diagonal;
	A.diag(-1) = -r_inner_diagonal;

	// Checkpoint
	std::cout << "Adding the r_inner_diagonal to A: " << "\n";
	std::cout << A	<< std::endl;

	// Defining two vectos with cx_vec for the sake of complex numbers.
	arma::vec a_vec(L, arma::fill::none);
	a_vec.fill(5);

	std::cout << "Diagonal vector a: " << "\n";
	std::cout << a_vec << std::endl;

	// Adding the a and b vectors to the A and B matrices
	A.diag(0) = a_vec;

	std::cout << "Matrix A " << "\n";
	A.print();
}
