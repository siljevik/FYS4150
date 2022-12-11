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
	int 	M = 5; 			// M-2 = 3 Here
	double 	h = 1/M; 		// Stepsize (since x_max and y_max will be 1 due to boudary conditions)
	double 	T = 100;			// Total time (T timeunits?? s??)
	double 	dt = 1;			// Timestep (T timeunits?? s??)
	double 	r_val = 1.0; 	// SANDER: wat the heck is this?
	int 	L = pow(M-2,2); // Length of A and B  (LxL)

	/*====================================================*/
    /*~~~~ Making matrices filled with ones and zeros ~~~~*/
    /*====================================================*/
	arma::mat U_n(M,M, arma::fill::zeros); 	// Matrix with elements u_ij (Space where particles move?)
	arma::mat V(M, M, arma::fill::ones);	// Matrix with elements v_ij (Potentials)
	arma::mat A(L, L, arma::fill::zeros);
	arma::mat B(L, L, arma::fill::zeros);

	// Define an empty vector that goes through a function that returns
	// a full vector dependent on indices (i,j) that we can run through with a loop
	//cout << "Vector V: \n" << V;
	//cout << head.index_translator(M) << endl;

	// TRUR ME BURDE HA NOKE INNI EIN TIMELOOP


	////////////// SILJE TESTYTESTS //////////////
	funcs.matrix_filler( M, r_val, L, A, B);
	funcs.diagonal_fill_AB(M, h, dt, L, V, A, B);
	arma::vec b = funcs.Bu_b(M, L, V, B);
	arma::vec u_n_one = funcs.Au_b(A,b);


 	return 0;
}





/*
     |\ | |  ||\ \ /(_~     |~)|_~|\/||_~|\/||~)|_~|~)
     |~\|_|/\||~\ | ,_)     |~\|__|  ||__|  ||_)|__|~\

        \ //~\| |    |\ |~)|_~    | ||\ ||/~\| ||_~
         | \_/\_/    |~\|~\|__    \_/| \||\_X\_/|__

      (J U S T   L I K E   E V E R Y O N E   E L S E)
      _____         _____         _____         _____
    .'     '.     .'     '.     .'     '.     .'     '.
   /  o   o  \   /  o   o  \   /  o   o  \   /  o   o  \
  |           | |           | |           | |           |
  |  \     /  | |  \     /  | |  \     /  | |  \     /  |
   \  '---'  /   \  '---'  /   \  '---'  /   \  '---'  /
jgs '._____.'     '._____.'     '._____.'     '._____.'
      _____         _____         _____         _____
    .'     '.     .'     '.     .'     '.     .'     '.
   /  o   o  \   /  o   o  \   /  o   o  \   /  o   o  \
  |           | |           | |           | |           |
  |  \     /  | |  \     /  | |  \     /  | |  \     /  |
   \  '---'  /   \  '---'  /   \  '---'  /   \  '---'  /
    '._____.'     '._____.'     '._____.'     '._____.'
      _____         _____         _____         _____
    .'     '.     .'     '.     .'     '.     .'     '.
   /  o   o  \   /  o   o  \   /  o   o  \   /  o   o  \
  |           | |           | |           | |           |
  |  \     /  | |  \     /  | |  \     /  | |  \     /  |
   \  '---'  /   \  '---'  /   \  '---'  /   \  '---'  /
    '._____.'     '._____.'     '._____.'     '._____.'
      _____         _____         _____         _____
    .'     '.     .'     '.     .'     '.     .'     '.
   /  o   o  \   /  o   o  \   /  o   o  \   /  o   o  \
  |           | |           | |           | |           |
  |  \     /  | |  \     /  | |  \     /  | |  \     /  |
   \  '---'  /   \  '---'  /   \  '---'  /   \  '---'  /
    '._____.'     '._____.'     '._____.'     '._____.'
      _____         _____         _____         _____
    .'     '.     .'     '.     .'     '.     .'     '.
   /  o   o  \   /  o   o  \   /  o   o  \   /  o   o  \
  |           | |           | |           | |           |
  |  \     /  | |  \     /  | |  \     /  | |  \     /  |
   \  '---'  /   \  '---'  /   \  '---'  /   \  '---'  /
    '._____.'     '._____.'     '._____.'     '._____.'
*/