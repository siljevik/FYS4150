#include <iostream>
#include <armadillo>

// problem 4
using namespace arma;
int prob4() {



    /////////////// FROM PROBLEM 2

    // big N to adjust the size of the NxN matrix A
        int N = 6;
    // minimum and maximum points
        double x_min = 0.0;
        double x_max = 1.0;
        int  n = N + 1;
        double dx = x_max - x_min;
        double h = dx*1./n;
    // create an empty matrix A, using Armadillo
        mat A(N, N, fill::zeros);
    // Creating two loops to fill the tridiagonal matrix with values:
    // values, a ,diag.
        for (int i=0; i<N; i++){
        A(i, i) =  2./(h*h);
        }
    // values,c ,diag.
        for (int i=0; i<(N-1) ; i++){
        A(i+1, i) = -1./(h*h);
        A(i, i+1) = -1./(h*h);
        }



    //############################################          Nødvendig??
    /////////////// FROM SECOND HINT IN PROBLEM 5 

    // Symmetrize the matrix by reflecting the upper triangle to lower triangle (just in case)
    A = arma::symmatu(A); 
    //############################################



    /////////////// FROM CODE SNIPPET AT END OF PROJECT DESCRIPTION

    // Determine the the max off-diagonal element of a symmetric matrix A
    // - Saves the matrix element indicies to k and l 
    // - Returns absolute value of A(k,l) as the function return value
    double max_offdiag_symmetric(const arma::mat& A, int& k, int& l);

    // Performs a single Jacobi rotation, to "rotate away"
    // the off-diagonal element at A(k,l).
    // - Assumes symmetric matrix, so we only consider k < l
    // - Modifies the input matrices A and R
    void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l);

    // Jacobi method eigensolver:
    // - Runs jacobo_rotate until max off-diagonal element < eps
    // - Writes the eigenvalues as entries in the vector "eigenvalues"
    // - Writes the eigenvectors as columns in the matrix "eigenvectors"
    //   (The returned eigenvalues and eigenvectors are sorted using arma::sort_index)
    // - Stops if it the number of iterations reaches "maxiter"
    // - Writes the number of iterations to the integer "iterations"
    // - Sets the bool reference "converged" to true if convergence was reached before hitting maxiter
    void jacobi_eigensolver(const arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors, 
                            const int maxiter, int& iterations, bool& converged);



    /////////////// FROM PROBLEM 2

    // printing out the tridiagonal matrix
        A.print(std::cout);
    // eigenvalues in order
        eigvalues.print(std::cout);
    // eigevectors in a matrix
        eigvectors.print(std::cout);

return 0;
}