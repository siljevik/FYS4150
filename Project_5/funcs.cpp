// =========================== //
// See bottom for a good laugh //
// =========================== //
#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <random>
#include <vector>
#include <map>

#include "funcs.hpp"

using namespace std;



/*===================================*/
/*~~~~~~     Vector Filler     ~~~~~~*/     // Problem 2.1
/*===================================*/
arma::cx_vec funcs::vector_filler(int M, arma::cx_mat V){
    // Defining an empty vector
    arma::cx_vec un_vec(pow(M-2,2), arma::fill::none);

	// fill the vector with vectors
    for(int i_p = 0; i_p < M-2; i_p++){
		int column = i_p;
		int c = column*(M-2);
		for(int j = 0; j < M-2; j++){
			// Append uij values into the uij_vec vector which will be appended
			// into un_vec
			int plc = c+j;
			un_vec(plc) = V(column,j);
        }//end of j-loop
    }// end of i-loop
	cout << "Vector un_vec in vector_filler: \n" << un_vec;
    return un_vec; // Returns the vector
}



/*======================================*/
/*~~~~~~     Index Translator     ~~~~~~*/     // Problem 2.1
/*======================================*/
void funcs::index_translator(int M, int k, int & i_p, int & j){
	// Test
	int column_from_k = (k+1)/(M-2); // int/int will make 5/2 = 2, or 6/7 = 0
	i_p = k-(column_from_k*(M-2))+1;
	j = column_from_k + 1;
}



/*===================================*/
/*~~~~~~     Matrix Filler     ~~~~~~*/      // Problem 2.2
/*===================================*/
void funcs::matrix_filler(int M, arma::cx_double r_val, int L, arma::cx_mat & A, arma::cx_mat & B){

	// Defining vector containing value r
	arma::cx_vec r_outer_diagonal( (M-2)*(M-3), arma::fill::none );// M-3 because we remove 1 position to shift the position in A/B.
	r_outer_diagonal.fill( r_val );

	// Checkpoint
	//std::cout << "Outer diagonal vector: " << "\n";
	//std::cout << r_outer_diagonal << std::endl;

	// Want to insert the diaginal matrices containing r into A and B
	A.diag(M-2) = -r_outer_diagonal;
	A.diag(2-M) = -r_outer_diagonal;
	B.diag(M-2) = r_outer_diagonal;
	B.diag(2-M) = r_outer_diagonal;

	// Checkpoint
	//std::cout << "Adding outer diagonal to B: " << "\n";
	//std::cout << B << std::endl;


	// Inner diagonal vectors containng r and for each 3rd position,0.
    arma::cx_vec r_inner_diagonal( (M-1)*(M-3), arma::fill::none );// M-3 because we remove 1 position to shift the position in A/B.
	r_inner_diagonal.fill(r_val);
	for(int k=0; k < (M-1)*(M-3); k++){
        if( (k+1)%3 == 0 )
        {
            r_inner_diagonal[k] = 0;
        } // end of if-statement
    }// end of k-loop

	// Checkpoint
	//std::cout << "r_inner_diagonal vector: " << "\n" ;
	//std::cout << r_inner_diagonal	<< std::endl;

	A.diag(1)  = -r_inner_diagonal;
	A.diag(-1) = -r_inner_diagonal;
    B.diag(1)  =  r_inner_diagonal;
    B.diag(-1) =  r_inner_diagonal;

	// Checkpoint
	//std::cout << "Adding the r_inner_diagonal to B: " << "\n";
	//std::cout << B	<< std::endl;

	// Defining two vectos with cx_vec for the sake of complex numbers.
	//arma::vec a_vec(L, arma::fill::none);
	//arma::vec b_vec(L, arma::fill::none);
	//a_vec.fill(5);
	//b_vec.fill(77);

	//std::cout << "Diagonal vector b: " << "\n";
	//std::cout << b_vec << std::endl;

	// Adding the a and b vectors to the A and B matrices
	//A.diag(0) = a_vec;
	//B.diag(0) = b_vec;

	//std::cout << "Matrix B " << "\n";
	//B.print();
}



/*================================================*/
/*~~~~~~~~   Diagonal filler of A and B   ~~~~~~~~*/     // Problem 2.3
/*================================================*/
// FROM PROBLEM 2: Now you are ready to write a function for your program 
// that, using inputs M, h,  and the matrix V as input, can fill two  
// matrices A and B and  according to the above pattern (point before this one)
void funcs::diagonal_fill_AB(int M, arma::cx_double h, arma::cx_double dt, int L, arma::cx_mat V,arma::cx_mat & A, arma::cx_mat & B){
	// Making the vector
	arma::cx_vec un_vec = vector_filler(M,V);
	// Calling i and j with the extension _plc to not confuse place i with
	// the complex number i.
	int i_plc;
	int j_plc;
	arma::cx_double icx = sqrt(-1);
	arma::cx_double r = (icx*dt)/(2.0*(pow(h,2)));

	for(int k = 0; k < L; k++)
	{
		//Finding the indices i and j
		index_translator(M, k, i_plc, j_plc); // don't use the & when using the function
		// vij is element place (i,j) in matrix V
		arma::cx_double vij = V(i_plc,j_plc);
		// Complex number i
		
		// Calculating a_k and b_k ----- SHOULD WE USE INT??? IDK
		arma::cx_double a_k = 1.0 + (4.0*r) + ((icx*dt)/2.0)*vij;
		arma::cx_double b_k = 1.0 - (4.0*r) - ((icx*dt)/2.0)*vij;
		
		A(k,k) = a_k;
		B(k,k) = b_k;
	}
	//cout << "\n Matrix A: \n" << A;
	cout << "\n Matrix B: \n" << B;
}



/*=====================================================*/
/*~~~~~~~~   Calculating b from the Bu^n = b   ~~~~~~~~*/      // Problem 3.1
/*=====================================================*/
arma::cx_vec funcs::Bu_b(int M, int L, arma::cx_mat V, arma::cx_mat B)
{
	// Random vector
	arma::cx_vec un_vec = vector_filler(M,V);
	//cout << "Vector un_vec in Bu_b: \n" << un_vec;//size(un_vec);
	//We need to transpose the vector (flip it 90 degrees)
	arma::cx_vec b = B*un_vec;
	cout << "\n\n VECTOR b: \n" << b;
	return b;
}



/*=======================================================*/
/*~~~~~~~~   Calculating the X from the AX = b   ~~~~~~~~*/     // Problem 3.2
/*=======================================================*/
arma::cx_vec funcs::Au_b(arma::cx_mat A, arma::cx_vec b)
{
	// Using Armadillo solve to solve for x in the 
	arma::cx_vec X = solve(A, b);
	cout << "\n Vector u^(n+1): \n" << X;
	return X;
}



/*===================================*/
/*~~~~~~~~   Initial state   ~~~~~~~~*/     // Problem 4
/*===================================*/
void funcs::initial_u(int M, double h, int L, arma::cx_vec U_0) // U_0 er initial matrix (hvor partiklene er????????????? halp)
{	// =========================================================================  //
	// Denne er ikke helt ferdig SOS
	// - Vet ikke hva x_c og y_c er (centre of initial wave packet???)
	// - Vet ikke hva sigma_x/y er (initial widths of the wave peacket???)
	// - Vet ikke hva p_x/y er heller (Wave packet momenta?? hææææææææææææææææææ)
	// ========================================================================= //

	// Declearing our variables here as doubles.
	arma::cx_double x_c, y_c; // coordinates of centre of initial wave packet
	arma::cx_double sigma_x, sigma_y; // Initial width of wave packet in x/y-direction
	arma::cx_double p_x, p_y; // Wave packet momenta

	// Since we will use the positions in calculations, we need y and x to be 
	// of the double-type variable
	for (arma::cx_double y = 0; y < M; y+=h){ // Since y_j = j*h;
		// To save time, let's do the y-calculations before next for-loop
		arma::cx_double division_y = ((y-y_c)^2)/(2*sigma_y^2);
		arma::cx_double ip_y = p_y*(y-y_c);//multiply with i (as in sqrt(-1))

		for (arma::cx_double x = 0; x < M; x+=h){ // Since x_i = i*h;
			arma::cx_double division_x = ((x-x_c)^2)/(2*sigma_x^2);
			arma::cx_double ip_x = p_x*(x-x_c);//multiply with i (as in sqrt(-1))

			// Unnormalized Gaussian wave packet -- Wave function then?
			U_0(x,y) = exp(-(division_x)-(division_y)+(ip_x)+(ip_y));
		} // End of x-loop
	} // End of y-loop

}




/*=======================================*/
/*~~~~~~~~   Initial potential   ~~~~~~~~*/     // Problem 5
/*=======================================*/

void funcs::initial_V()
{
	// IDK WHAT I'M DOING HALP

	// Code here
	arma::cx_double wt = 0.02; 		// Wall thickness
	arma::cx_double x_c = 0.5; 		// x-centre position
	arma::cx_double y_c = 0.5;		// y-centre
	arma::cx_double wts_y = 0.05; 	// Thickness of wall piece separating the two slits (y-distance between the inner edges of the two slits)
	arma::cx_double so_y = 0.05; 	// Slit opening (y-direction)
	// Ensure that the slit setup is symmetric around y_c (wall between slits in the middle of y_c)

}


// Fight Bugs                      |     |
//                                 \\_V_//
//                                 \/=|=\/
//                                  [=v=]
//                                __\___/_____
//                               /..[  _____  ]
//                              /_  [ [  M /] ]
//                             /../.[ [ M /@] ]
//                            <-->[_[ [M /@/] ]
//                           /../ [.[ [ /@/ ] ]
//      _________________]\ /__/  [_[ [/@/ C] ]
//     <_________________>>0---]  [=\ \@/ C / /
//        ___      ___   ]/000o   /__\ \ C / /
//           \    /              /....\ \_/ /
//        ....\||/....           [___/=\___/
//       .    .  .    .          [...] [...]
//      .      ..      .         [___/ \___]
//      .    0 .. 0    .         <---> <--->
//   /\/\.    .  .    ./\/\      [..]   [..]
//  / / / .../|  |\... \ \ \    _[__]   [__]_
// / / /       \/       \ \ \  [____>   <____]