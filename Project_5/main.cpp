///////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                  		 //
// Compile with: g++ -O2 main.cpp -std=c++11 -I include -o main.exe -lomp -larmadillo -lsuperlu    //
//                                                                                  		 //
// Run with: ./main.exe                                                              	 	 //
//                                                                                  		 //
///////////////////////////////////////////////////////////////////////////////////////////////
// SEE BOTTOM FOR CHRISTMAS SURPRISE <3 //
//////////////////////////////////////////
#include <iostream>
#include <armadillo>
#include <complex>
#include <vector>

#include "funcs.hpp"
#include "funcs.cpp"

using namespace std;

int main(){

	funcs funcs;

	/*====================================*/
    /*~~~~ Constants, Variables, etc. ~~~~*/
    /*====================================*/
	// === From problem 5: === //
	arma::cx_double wtx = 0.02; 	// Wall thickness x-dir
	arma::cx_double x_c = 0.5; 		// x-centre position
	arma::cx_double y_c = 0.5;		// y-centre
	arma::cx_double wts_y = 0.05; 	// Thickness of wall piece separating the two slits (y-distance between the inner edges of the two slits)
	arma::cx_double so_y = 0.05; 	// Slit opening (y-direction)
	// === From problem 7: === //
	arma::cx_double h_cx = 0.05; //////////////////// ENDRE TIL RIKTIG ETTER TESTING
	arma::cx_double dt = 2.5*pow(10,-5);
	arma::cx_double T = 0.008;
	arma::cx_double sigma_x = 0.05;
	arma::cx_double sigma_y = 0.05;
	arma::cx_double p_x = 200;
	arma::cx_double p_y = 0;
	arma::cx_double v_0 = 0;
	// === From problem 8 === //
	/*
	arma::cx_double h_cx = 0.005;
	arma::cx_double dt = 2.5*pow(10,-5);
	arma::cx_double T = 0.002;
	arma::cx_double sigma_x = 0.05;
	arma::cx_double sigma_y = 0.05;
	arma::cx_double p_x = 200;
	arma::cx_double p_y = 0;
	arma::cx_double v_0 = 0;
	*/ 
	// Stuff we find using the information above
	double h = h_cx.real();
	arma::cx_double M_cx = 1/h; // Test?
	int M = M_cx.real();
	int L = pow(M-2,2); // Length of A and B 
	// Testing for r = (i*dt)/(2.0*(h^2))
	arma::cx_double icx(0.0,1.0);
	arma::cx_double r_val= icx*(dt)/(2.0*(h*h));



	/*====================================================*/
    /*~~~~ Making matrices filled with ones and zeros ~~~~*/
    /*====================================================*/
	arma::cx_mat U_n(M,M, arma::fill::zeros); 	// Matrix with elements u_ij (Space where particles move)
	arma::sp_cx_mat spU(U_n);					// Making the matrix sparse to save time
	arma::cx_mat V(M, M, arma::fill::zeros);		// Matrix with elements v_ij (Potentials)
	arma::sp_cx_mat spV(V);
	arma::cx_mat A(L, L, arma::fill::zeros);
	arma::sp_cx_mat spA(A);
	arma::cx_mat B(L, L, arma::fill::zeros);
	arma::sp_cx_mat spB(B);


	
	/*====================================================*/
    /*~~~~ Making matrices filled with ones and zeros ~~~~*/     // Problem 6
    /*====================================================*/
// 1. Simulation parameters given above

// 2. Setting up the potential matrix V:
	funcs.double_slit(h, v_0, M, spV);

// 3. Setting up the initial state matrix U^0:
	arma::cx_vec u0 = funcs.vector_filler(M, spU);
	funcs.initial_u(M,h,L,u0, spU, x_c, y_c, sigma_x, sigma_y, p_x, p_y);
	
// 4. Setting up the matrices A and B:
	funcs.matrix_filler(M, r_val, L, spA, spB);			// Inserting the r
	funcs.diagonal_fill_AB(M, h, dt, L, spV, spA, spB);   // Diagonals (a_k and b_k)
	arma::cx_vec b = funcs.Bu_b(M, L, spU, spB);
	
	arma::cx_vec u_np1 = funcs.Au_b(spA,b); // Den skjønner ikkje 'spsolve()' :(

// 5. Running a loop over time steps to store each new state U^n	
	double time_step = (double)dt.real();//(double)dt.real();
	double total_time = (double)T.real();
	int total_time_steps = (total_time/time_step)/10;
	cout << "total time steps: " << total_time_steps;

	// Creating and opening a .txt-file
	ofstream datafile;
	datafile.open("cx_cube.txt", ofstream::out | ofstream::trunc);

	//Creating the cx_cube
	arma::cx_cube cubeboi(M,M,total_time_steps+1);
	cubeboi.slice(0) = spU; // Initial matrix inserted to the cube
	
	// Then we are ready to run the time-loop:
	for (int t = 1; t <= total_time_steps; t++)
	{
		arma::cx_vec b = funcs.Bu_b(M, L, spU, spB);
		arma::cx_vec u_t = funcs.Au_b(spA,b);
		// Not initial U, but U^n
		funcs.initial_u(M,h,L, u_t, spU, x_c, y_c, sigma_x, sigma_y, p_x, p_y);
		cubeboi.slice(t) = spU;
	}
	// Adding cube to txt file:
	datafile << cubeboi;
	datafile.close(); // Closing datafile
	

 	return 0;
}



//        SEASONS GREETINGS and the best for the NEW YEAR!
//                                                ._...._.
//                        \ .*. /                .::o:::::.
//                         (\o/)                .:::'''':o:.
//                          >*<                 `:}_>()<_{:'
//                         >0<@<             @    `'//\\'`    @
//                        >>>@<<*          @ #     //  \\     # @
//                       >@>*<0<<<      .__#_#____/'____'\____#_#_.
//                      >*>>@<<<@<<     [_________________________]
//                     >@>>0<<<*<<@<     |=_- .-/\ /\ /\ /\--. =_-|
//                    >*>>0<<@<<<@<<<    |-_= | \ \ \ \ \ \\-|-_=-|
//                   >@>>*<<@<>*<<0<*<   |_=-=| / // // // / |_=-_|
//     \*/         //0>>*<<@<>0><<*<@<<  |=_- | `-'`-'`-'`-' |=_=-|
// .___\U//__.    >*>>@><0<<*>>@><*<0<<  | =_-| o          o |_==_|
//  \ | | \  |  >@>>0<*<<0>>@<<0<<<*<@<  |=_- | !     (    ! |=-_=|
//   \| | _(UU)_ >((*))_>0><*<0><@<<<0<*<|-,-=| !    ).    ! |-_-=|
//  \ \| || / //||.*.*.*.|>>@<<*<<@>><0<<|=_,=| ! __(:')__ ! |=_==|
//  \_|_|&&_// ||*.*.*.*|_\db//_   (\_/)-|     /^\=^=^^=^=/^\| _=_|
//   ""|'.'.'.|~~|.*.*.*|      |  =('Y')=|=_,//.------------.\\_,_|
//     |'.'.'.|  |^^^^^^|______|  ( ~~~ )|_,_/(((((((())))))))\_,_|
//     ~~~~~~~ ""       `------'  `w---w'|_____`------------'_____|
// ________________________________________________________________