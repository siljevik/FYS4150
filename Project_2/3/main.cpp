#include <iostream>
#include <armadillo>
#include <assert.h>

// problem 3
using namespace arma;
// declaration
mat create_symmetric(int n, double a, double d);
double max_offdiag_symmetric(const mat&A, int& k, int&l);



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
	cout << max_offdiag_symmetric( create_symmetric(n,a,d) ,k,l) << endl;
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
// - The matrix indices of the max element are returned by writing to the
//   int references k and l (row and column, respectively)
// - The value of the max element A(k,l) is returned as the function
//   return value


double max_offdiag_symmetric(const mat& A, int& k, int& l){
  // Getting the size of matrix A
	int N;
	A.n_rows;

  // Checking the consistancy; if the matrix is square and larger than 1x1, with armadillo. Prints out "Is Matrix A square?: true" if true
//	cout << "Is Matrix A square?:" << endl;
//	cout << A.is_square() << endl;

  // The standard function 'assert' from <assert.h> can be useful for quick checks like this
  // during the code development phase. Use it like this: assert(some condition),
  // e.g assert(a==b). If the condition evaluates to false, the program is killed with
  // an assertion error. More info: https://www.cplusplus.com/reference/cassert/assert/
	// quick check if the matrix
	assert (A.is_square() == true);
  // Initialize references k and l to the first off-diagonal element of A
  // Initialize a double variable 'maxval' to A(k,l). We'll use this variable
  // to keep track of the largest off-diag element.
//	double maxval = A.max();

  // Loop through all elements in the upper triangle of A (not including the diagonal)
  // When encountering a matrix element with larger absolute value than the current value of maxval,
  // update k, l and max accordingly.

// looping through the upper triangle of the NxN matrix A.
// When encountering a matrix element with a larger absolute value than the current value of maxval,
// k, l and maxval will be updated

	double max;
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
            cout << max;
}

