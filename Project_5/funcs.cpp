// ========================== //
// See bottom for a good time //
// ========================== //
#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <complex>
#include <random>
#include <vector>
#include <map>

#include "funcs.hpp"

using namespace std;



/*===================================*/
/*~~~~~~     Vector Filler     ~~~~~~*/     // Problem 2.1
/*===================================*/
arma::cx_vec funcs::vector_filler(int M, arma::sp_cx_mat V){
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
	//cout << "Vector un_vec in vector_filler: \n" << un_vec;
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
void funcs::matrix_filler(int M, arma::cx_double r_val, int L, arma::sp_cx_mat & A, arma::sp_cx_mat & B){

	// Defining vector containing value r
	arma::cx_vec r_outer_diagonal( (M-2)*(M-3), arma::fill::none );// M-3 because we remove 1 position to shift the position in A/B.
	r_outer_diagonal.fill( r_val );

	// Want to insert the diaginal matrices containing r into A and B
	A.diag(M-2) = -r_outer_diagonal;
	A.diag(2-M) = -r_outer_diagonal;
	B.diag(M-2) = r_outer_diagonal;
	B.diag(2-M) = r_outer_diagonal;

	// Inner diagonal vectors containng r and for each 3rd position,0.
    arma::cx_vec r_inner_diagonal( (M-1)*(M-3), arma::fill::none );// M-3 because we remove 1 position to shift the position in A/B.
	r_inner_diagonal.fill(r_val);
	for(int k=0; k < (M-1)*(M-3); k++){
        if( (k+1)%3 == 0 )
        {
            r_inner_diagonal[k] = 0;
        } // end of if-statement
    }// end of k-loop

	A.diag(1)  = -r_inner_diagonal;
	A.diag(-1) = -r_inner_diagonal;
    B.diag(1)  =  r_inner_diagonal;
   	B.diag(-1) =  r_inner_diagonal;
}



/*================================================*/
/*~~~~~~~~   Diagonal filler of A and B   ~~~~~~~~*/     // Problem 2.3
/*================================================*/
// FROM PROBLEM 2: Now you are ready to write a function for your program 
// that, using inputs M, h,  and the matrix V as input, can fill two  
// matrices A and B and  according to the above pattern (point before this one)
void funcs::diagonal_fill_AB(int M, arma::cx_double h, arma::cx_double dt, int L, arma::sp_cx_mat V,arma::sp_cx_mat & A, arma::sp_cx_mat & B){
	
	// Making the vector
	arma::cx_vec un_vec = vector_filler(M,V);
	
	// Calling i and j with the extension _plc to not confuse place i with
	// the complex number i.
	int i_plc;
	int j_plc;
	arma::cx_double icx(0.0,1.0); // Complex i
	arma::cx_double r = (icx*dt)/(2.0*(pow(h,2)));

	for(int k = 0; k < L; k++)
	{
		//Finding the indices i and j
		index_translator(M, k, i_plc, j_plc); // don't use the & when using the function
		
		// vij is element place (i,j) in matrix V
		arma::cx_double vij = V(i_plc,j_plc);

		// Calculating a_k and b_k
		arma::cx_double a_k = 1.0 + (4.0*r) + ((icx*dt)/2.0)*vij;
		arma::cx_double b_k = 1.0 - (4.0*r) - ((icx*dt)/2.0)*vij;

		// Inserting them into the diagonal spots
		A(k,k) = a_k;
		B(k,k) = b_k;
	}
}



/*=====================================================*/
/*~~~~~~~~   Calculating b from the Bu^n = b   ~~~~~~~~*/      // Problem 3.1
/*=====================================================*/
arma::cx_vec funcs::Bu_b(int M, int L, arma::sp_cx_mat U, arma::sp_cx_mat B)
{
	// First we make the vector containing all elements from U
	arma::cx_vec un_vec = vector_filler(M,U);
	// Then we can use it to solve for vector b
	arma::cx_vec b = B*un_vec;
	return b;
}



/*=======================================================*/
/*~~~~~~~~   Calculating the X from the AX = b   ~~~~~~~~*/     // Problem 3.2
/*=======================================================*/
arma::cx_vec funcs::Au_b( arma::sp_cx_mat A, arma::cx_vec b)
{
	// Using Armadillo solve to solve for x in the Ax=b equation, where x here is u_np1
	arma::cx_vec u_np1 = spsolve(A, b);
	return u_np1;
}



/*===================================*/
/*~~~~~~~~   Initial state   ~~~~~~~~*/     // Problem 4
/*===================================*/
void funcs::initial_u(int M, double h, int L, arma::cx_vec u0, arma::sp_cx_mat & U_n, arma::cx_double x_c, arma::cx_double y_c,
		arma::cx_double sigma_x, arma::cx_double sigma_y, arma::cx_double p_x, arma::cx_double p_y)
{
	// Since we want to use M in our for-loops, we convert it to a double.
	// Adn we will use the positions in calculations, we need y and x to be 
	// of the double-type variable
	double M_d = (double)M;

	// Defining a normalized vector U_n
	arma::sp_cx_mat norm_U_n;
	// To save time, we will do some calcuations outside the for-loop when we can:
	arma::cx_double under_y_division = 2.0*pow(sigma_y, 2);
	arma::cx_double under_x_division = 2.0*pow(sigma_x, 2);
	for (double y = 0; y < M_d; y+=h){ // Since y_j = j*h;
		// To save time, let's do the y-calculations before next for-loop
		arma::cx_double division_y = (pow((y-y_c),2))/under_y_division;
		arma::cx_double ip_y = p_y*(y-y_c);//multiply with i (as in sqrt(-1))

		for (double x = 0; x < M_d; x+=h){ // Since x_i = i*h;
			arma::cx_double division_x = (pow((x-x_c),2))/under_x_division;
			arma::cx_double ip_x = p_x*(x-x_c);//multiply with i (as in sqrt(-1))

			// Unnormalized Gaussian wave packet -- Wave function then?
			U_n(x,y) = exp(-(division_x)-(division_y)+(ip_x)+(ip_y));

			arma::sp_cx_mat conj_Un = arma::conj(U_n);
			// Normalized Gaussian
			norm_U_n = conj_Un * U_n(x,y);
		} // End of x-loop
	} // End of y-loop

}




/*=======================================*/
/*~~~~~~~~   Initial potential   ~~~~~~~~*/     // Problem 5
/*=======================================*/
void funcs::double_slit(double h, arma::cx_double v0, int M, arma::sp_cx_mat & V)
{
	// We need to set all points that is not wall or on the borders to v0.
	// So, here we just ignore/don't do anything to the walls:D
	//double h = h_cx.real(); // First to do stuff with h in the loops, we convert it to 

	// Area before slit (left side)
	int stopx_1 = 0.49/h; // Doing these outside for-loops to save time
	for(int i = 1; i < M-1; i++){			// y
		for(double j = 1; j < stopx_1; j+=1){ // x
			V(i,j) = v0;
		}
	}

	// Area in top slit
	int starty_2 = 0.425/h;
	int stopy_2  = (0.475+h)/h; 
	int startx_2 = 0.49/h;
	int stopx_2  = (0.51+h)/h;
	for(int k = starty_2; k < stopy_2-1; k++){  // kkk ILLUMINATI
		for(int l = startx_2; l < stopx_2; l++){
			V(k,l) = v0;
		}
	}

	// Area in bottom slit
	int starty_3 = 0.525/h;
	int stopy_3  = (0.575+h)/h;
	int startx_3 = 0.49/h;
	int stopx_3  = (0.51+h)/h;
	for(int m = starty_3; m < stopy_3-1; m++){
		for(int n = startx_3; n < stopx_3; n++){
			V(m,n) = v0;
		}
	}

	// Area after slit (right side)
	int startx_4 = (0.51+h)/h;
	int stopx_4  = ((1-h)/h)+1;
	for(int o = 1; o < M-1; o++){
		for(int p = startx_4; p < stopx_4; p++){
			V(o,p) = v0;
		}
	}
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
