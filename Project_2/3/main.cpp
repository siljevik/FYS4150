#include <iostream>
#include <armadillo>
#include <assert.h>

// problem 3
using namespace arma;
// declaration
mat create_symmetric(int n, double a, double d);
double max_offdiag_symmetric(const mat&A, int& k, int&l, int& N);



//----- main function -----
int main(){
// big N to adjust the size of the NxN matrix A
	int N = 4;
// minimum and maximum points
	double x_min = 0.0;
	double x_max = 1.0;

	int  n = N + 1;
	double dx = x_max - x_min;
	double h = dx*1./n;

// tridiagonal values
	double a = -1./(h*h);
	double d =  2./(h*h);

// k and l
	int k = 0;
	int l = N;

// print matrix
	cout << "Our Matrix A:" << endl;
	cout << create_symmetric(n,a,d) << endl;
	cout << "With the largest off-diagonal element(in abs.):";
	cout << max_offdiag_symmetric( create_symmetric(n,a,d) ,k,l, N) << endl;
return 0;
}


mat create_symmetric(int n, double a, double d){
	int N = n-1;
// create an empty matrix A, using Armadillo
	mat A(N, N, fill::randn);

// Creating two loops to fill the tridiagonal matrix with values:
// values, d ,diag.
	for (int i=0; i<N; i++)
		{
 		A(i, i) =  d;
		}

// values,a ,diag.
	for (int i=0; i<(N-1) ; i++)
		{
 		A(i+1, i) = a;
 		A(i, i+1) = a;
		}

	return A;
}


// A function that finds the max off-diag element of a symmetric matrix A.
double max_offdiag_symmetric(const mat& A, int& k, int& l, int& N){
	// quick check if the matrix is symmetric with the .is_square() function from armadillo
	assert (A.is_square() == true);
	
	double max;	// maximum value of A(k,l)
            for (int i = N-1; i > -1; --i) // Looking at position A(0,N-1) to A(N-1,0)
            {
                for ( int j = i; j > -1; --j)
                {
                    double aij = fabs(A(i,j));
                    if ( aij > max)
                    {
                        max = aij; k = i; l = j;
                    }
                }
            }
            return max;
}

