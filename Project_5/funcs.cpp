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
	arma::mat B(L,L, arma::fill::zeros);

	// Defining vector containing value r
	arma::vec r_outer_diagonal( (M-2)*(M-3), arma::fill::none );// M-3 because we remove 1 position to shift the position in A/B.
	r_outer_diagonal.fill( r_val );

	// Checkpoint
	std::cout << "Outer diagonal vector: " << "\n";
	std::cout << r_outer_diagonal << std::endl;

	// want to insert the diaginal matrices containing r into A and B
	A.diag(M-2) = -r_outer_diagonal;
	A.diag(2-M) = -r_outer_diagonal;
	B.diag(M-2) = r_outer_diagonal;
	B.diag(2-M) = r_outer_diagonal;

	// Checkpoint
	std::cout << "Adding outer diagonal to B: " << "\n";
	std::cout << B << std::endl;


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
        B.diag(1)  = r_inner_diagonal;
        B.diag(-1) = r_inner_diagonal;

	// Checkpoint
	std::cout << "Adding the r_inner_diagonal to B: " << "\n";
	std::cout << B	<< std::endl;

	// Defining two vectos with cx_vec for the sake of complex numbers.
	arma::vec a_vec(L, arma::fill::none);
	arma::vec b_vec(L, arma::fill::none);
	a_vec.fill(5);
	b_vec.fill(77);

	std::cout << "Diagonal vector b: " << "\n";
	std::cout << b_vec << std::endl;

	// Adding the a and b vectors to the A and B matrices
	A.diag(0) = a_vec;
	B.diag(0) = b_vec;

	std::cout << "Matrix B " << "\n";
	B.print();
}

/*================================================*/
/*~~~~~~~~   Diagonal filler of A and B   ~~~~~~~~*/
/*================================================*/
// FROM PROBLEM 2: Now you are ready to write a function for your program 
// that, using inputs M, h,  and the matrix V as input, can fill two  
// matrices A and B and  according to the above pattern (point before this one)

void Header::diagonal_fill_AB(int M, int h, arma::mat V,arma::mat & A, arma::mat & B){
	// Making the vector
	arma::vec un_vec = vector_filler(V);
	// Calling i and j with the extension _plc to not confuse place i with
	// the complex number i.
	int i_plc;
	int j_plc;

	////////////////////////////////////////////////
	// For testing, matrisen lages i problem 5
	arma::mat V(M,M, arma::fill::ones);
	cout << "\n Matrix V: \n" << V ;
	//v_ij = V(i,j) -- V er matrise, vij er element i matrise
	////////////////////////////////////////////////

	int length_ks = pow((M-2),2)
	for(int k = 0; k < length_ks; k++)
	{
		//Finding the indices i and j
		index_translator(V, k, & i_plc, & j_plc);
		// vij is element place (i,j) in matrix V
		int vij = V(i_plc,j_plc);
		// Complex number i
		arma::cx_double i_comp = 1i;
		// Calculating a_k and b_k ----- SHOULD WE USE INT??? IDK
		int a_k = 1 + (4*r) + (i_comp*trekant_t/2)*vij;
		int b_k = 1 - (4*r) - (i_comp*trekant_t/2)*vij;

		A(k,k) = a_k;
		B(k,k) = b_k;
	}
	cout << "\n Matrix A: \n" << A;
	cout << "\n Matrix B: \n" << B;
}
