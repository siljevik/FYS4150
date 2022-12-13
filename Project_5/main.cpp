//////////////////////////////////////////////////////////////////////////////////////
//                                                                                  //
// Compile with: g++ -O2 main.cpp -std=c++11 -I include -o main.exe -larmadillo     //
//                                                                                  //
// Run with: ./main.exe                                                             //
//                                                                                  //
//////////////////////////////////////////////////////////////////////////////////////
// SEE BOTTOM FOR CHRISTMAS SURPRISE <3 //
//////////////////////////////////////////
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
	//int 	M = 5; 			// M-2 = 3 Here (LxL)
	//arma::cx_double 	h = 1/M; 		// Stepsize (since x_max and y_max will be 1 due to boudary conditions)
	//arma::cx_double 	T = 100;			// Total time (T timeunits?? s??)
	//arma::cx_double 	dt = 1;			// Timestep (T timeunits?? s??)
	// who dis
	 	// SANDER: wat the heck is this?
	
	////////////////////////////////////////////////////////////
	// Values that sould've been in a txt document or something
	// Will put into txt document if time
	// === From problem 5: === //
	arma::cx_double wtx = 0.02; 	// Wall thickness x-dir
	arma::cx_double x_c = 0.5; 		// x-centre position
	arma::cx_double y_c = 0.5;		// y-centre
	arma::cx_double wts_y = 0.05; 	// Thickness of wall piece separating the two slits (y-distance between the inner edges of the two slits)
	arma::cx_double so_y = 0.05; 	// Slit opening (y-direction)
	// === From problem 7: === //
	arma::cx_double h_cx = 0.005;
	arma::cx_double dt = 2.5*pow(10,-5);
	arma::cx_double T = 0.008;
	arma::cx_double sigma_x = 0.05;
	arma::cx_double sigma_y = 0.05;
	arma::cx_double p_x = 200;
	arma::cx_double p_y = 0;
	arma::cx_double v_0 = 0;
	

	// Stuff we find using the information above
	double h = h_cx.real();
	arma::cx_double M_cx = 1/h; // Test?
	int M = M_cx.real();
	int L = pow(M-2,2); // Length of A and B 
	// Testing for r = (i*dt)/(2.0*(h^2))
	arma::cx_double 	r_val= (sqrt(-1))*(dt)/(2.0*(h*h));

	/*====================================================*/
    /*~~~~ Making matrices filled with ones and zeros ~~~~*/
    /*====================================================*/
	arma::cx_mat U_n(M,M, arma::fill::zeros); 	// Matrix with elements u_ij (Space where particles move)
	arma::cx_mat V(M, M, arma::fill::none);		// Matrix with elements v_ij (Potentials)
	arma::cx_mat A(L, L, arma::fill::zeros);
	arma::cx_mat B(L, L, arma::fill::zeros);

	// Define an empty vector that goes through a function that returns
	// a full vector dependent on indices (i,j) that we can run through with a loop
	//cout << "Vector V: \n" << V;
	//cout << head.index_translator(M) << endl;

	// TRUR ME BURDE HA NOKE INNI EIN TIMELOOP


	////////////// SILJE TESTYTESTS //////////////
	//funcs.matrix_filler( M, r_val, L, A, B);
	//funcs.diagonal_fill_AB(M, h_cx, dt, L, V, A, B);
	//arma::cx_vec b = funcs.Bu_b(M, L, V, B);
	//arma::cx_vec u_n_one = funcs.Au_b(A,b);
	
	/*====================================================*/
    /*~~~~ Making matrices filled with ones and zeros ~~~~*/     // Problem 6
    /*====================================================*/
	// 1. Simulation parameters given above
	// 2. Setting up the potential matrix V:
	funcs.double_slit(h_cx, v_0, V);
	// 3. Setting up the initial state matrix U^0:
	funcs.initial_u(M,h,L,U_0, x_c, y_c, sigma_x, sigma_y, p_x, p_y);
	// 4. Setting up the matrices A and B:
	funcs.diagonal_fill_AB(M, h, dt, L, V, A, B);
	// 5. Running a loop over time steps to store each new state U^n	
	double time_step = 0.001;//(double)dt.real();
	double total_time = (double)T.real();
	for (double t = 0; t <= total_time; t+=time_step)
	{
		cout << "Hello World";
	}



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