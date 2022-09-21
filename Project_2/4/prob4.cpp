///////////////////////////////////////////
//  FOR SILJES EYES ONLY:                //
//  Compile with                         //
//  g++ -std=c++11 prob4.cpp -larmadillo //
//  Run with                             //
//  ./a.out                              //
///////////////////////////////////////////

#include <iostream>
#include <armadillo>


// problem 4
using namespace arma;
// Declarations
double offdiag(mat A, int p, int q, int N);
mat jacobi_rotate(mat A, mat R, int k, int l, int N);
//void jacobi_eigensolver(const arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors, 
//                            const int maxiter, int& iterations, bool& converged);


int main() {
    /////////////// FROM PROBLEM 2

    // big N to adjust the size of the NxN matrix A
        int N = 6;  // From Praoblem 4 b)
    // minimum and maximum points
        double x_min = 0.0;
        double x_max = 1.0;
        int  n = N + 1;
        double dx = x_max - x_min;
        double h = dx*1./n;
    // create an empty matrix A, using Armadillo
        mat A(N, N, fill::zeros);
        mat R(N, N, fill::zeros);  // Matrix R to be filled with the eigenvectors
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
    std::cout << "Matrix A: \n";
    std::cout << A;

    /////////////
    // we want to solve Av = lambda v using arma::eig_sym
	vec eigval; // vector containing eigenvalues
	mat eigvec; // eigenvector
	eig_sym(eigval, eigvec, A);
    // eigenvalues in order
    std::cout << "Eigen-values: \n" << eigval;
    // eigevectors in a matrix
    std::cout << "Eigen-vectors: \n"<< eigvec;
    ////////////


    int p = 0;
    int q = 0;
    int l = 0;
    int k = 0;
    double max;
    max = offdiag( A, p, q,N);
    std::cout << "The maximum (absolute) value of A is " << max << '\n';
    R = jacobi_rotate( A, R, k, l, N);
    //A.print(std::cout);
    std::cout << "The Jacobi rotated matrix is: \n";
    R.print(std::cout);
    return 0;
}



    ////////////// FROM CODE SNIPPET AT END OF PROJECT DESCRIPTION

        // Determine the the max off-diagonal element of a symmetric matrix A
        // - Saves the matrix element indicies to k and l 
        // - Returns absolute value of A(k,l) as the function return value
    double offdiag(mat A, int p, int q, int N) ///////////// OPPGAVE 3!!!
    {
            double max;
            for (int i = N-1; i > -1; --i) // Looking at position A(0,N-1) to A(N-1,0)
            {
                for ( int j = i; j > -1; --j)
                {
                    double aij = fabs(A(i,j));
                    if ( aij > max)
                    {
                        max = aij; p = i; q = j;
                    }
                }
            }
            return max; 
        }

//////////////////////////  T E S T  //////////////////////////
mat jacobi_rotate(mat A, mat R, int row, int col, int N)
    {
        int n = N - 1;
        mat new_A(N, N, fill::zeros); 
        for ( int row = 0; row < N; row++ ) 
        {
            for (int col = 0; col < N; col++)
            {
                if (row != col)
                {
                    int new_row = fabs(row - n);
                    int new_col = fabs(col - n);
                    new_A(new_row,new_col) = A(row,col);
                }
                else 
                {
                    new_A(row,col) = A(row,col);
                }
            }
        }
        return new_A;
    }


        /*double s, c;
        if ( A(k,l) != 0.0 ) {
            double t, tau;
            tau = (A(l,l) - A(k,k))/(2*A(k,l));
            if ( tau >= 0 ) {
                t = 1.0/(tau + sqrt(1.0 + tau*tau));
            } else {
                t = -1.0/(-tau + sqrt(1.0 + tau*tau));
            }
            c = 1/sqrt(1+t*t);
            s = c*t;
        } else {
            c = 1.0;
            s = 0.0;
        }
        double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
        a_kk = A(k,k);
        a_ll = A(l,l);
        A(k,k) = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll;
        A(l,l) = s*s*a_kk + 2.0*c*s*A(k,l) + c*c*a_ll;
        A(k,l) = 0.0; // hard-coding non-diagonal elements by hand
        A(l,k) = 0.0; // same here
        for ( int i = 0; i < N-1; i++ ) {
            if ( i != k && i != l ) {
                a_ik = A(i,k);
                a_il = A(i,l);
                A(i,k) = c*a_ik - s*a_il;
                A(k,i) = A(i,k);
                A(i,l) = c*a_il + s*a_ik;
                A(l,i) = A(i,l);
            
            // And finally the new eigenvectors
            r_ik = R(i,k);
            r_il = R(i,l);
            R(i,k) = c*r_ik - s*r_il;
            R(i,l) = c*r_il + s*r_ik;
        }
        }
        A.print(std::cout);
        R.print(std::cout);
        return R;
    }*/
//////////////////////////  T E S T  //////////////////////////





    // Performs a single Jacobi rotation, to "rotate away"
    // the off-diagonal element at A(k,l).
    // - Assumes symmetric matrix, so we only consider k < l
    // - Modifies the input matrices A and R
    /*mat jacobi_rotate(mat A, mat R, int k, int l, int N)
    {
        //R = mat(N, N, fill::zeros);
        
        double s, c;
        if ( A(k,l) != 0.0 ) 
        {    
            double t, tau;
            tau = (A(l,l) - A(k,k))/(2*A(k,l));

            if ( tau >= 0 )
                {t = 1.0/(tau + sqrt(1.0 + tau*tau));} 
            else 
                {t = -1.0/(-tau + sqrt(1.0 + tau*tau));}
            
            c = 1/sqrt(1+t*t);
            s = c*t;
        } 

        else 
        {
            c = 1.0;
            s = 0.0;
        }

        double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
        
        a_kk = A(k,k);
        a_ll = A(l,l);
        
        A(k,k) = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll;
        A(l,l) = s*s*a_kk + 2.0*c*s*A(k,l) + c*c*a_ll;
        
        A(k,l) = 0.0; // hard-coding non-diagonal elements by hand
        A(l,k) = 0.0; // same here
        
        for ( int i = 0; i < N-1; i++ ) 
        {
            if ( i != k && i != l ) 
            {
                a_ik = A(i,k);
                a_il = A(i,l);
                A(i,k) = c*a_ik - s*a_il;
                A(k,i) = A(i,k);
                A(i,l) = c*a_il + s*a_ik;
                A(l,i) = A(i,l);
                // And finally the new eigenvectors
                r_ik = A(i,k);
                r_il = A(i,l);
                R(i,k) = c*r_ik - s*r_il;
                R(i,l) = c*r_il + s*r_ik;
            }
        }
        
        A.print(std::cout);
        R.print(std::cout);
        return R;

    } */// end of function jacobi_rotate
 ///////////////// TEST
 

    // Jacobi method eigensolver:
    // - Runs jacobo_rotate until max off-diagonal element < eps
    // - Writes the eigenvalues as entries in the vector "eigenvalues"
    // - Writes the eigenvectors as columns in the matrix "eigenvectors"
    //   (The returned eigenvalues and eigenvectors are sorted using arma::sort_index)
    // - Stops if it the number of iterations reaches "maxiter"
    // - Writes the number of iterations to the integer "iterations"
    // - Sets the bool reference "converged" to true if convergence was reached before hitting maxiter
//    void jacobi_eigensolver(const arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors, 
//                            const int maxiter, int& iterations, bool& converged)



    /////////////// FROM PROBLEM 2

    // printing out the tridiagonal matrix
    //    A.print(std::cout);
    //    R.print(std::cout);
    // eigenvalues in order
        //eigvalues.print(std::cout);
    // eigevectors in a matrix
        //eigvectors.print(std::cout);

