///////////////////////////////////////////
//  FOR SILJES EYES ONLY:                //
//  Compile with                         //
//  g++ -std=c++11 prob4.cpp -larmadillo //
//  Run with                             //
//  ./a.out                              //
///////////////////////////////////////////

#include <iostream>
#include <armadillo>
using namespace arma;

// Declarations
double offdiag(mat& A, int& p, int& q, int N);
void jacobi_rotate(mat& A, mat& R, int k, int l, int N);
mat eigen_solver_loop(double tolerance, mat& A, mat& R, int N, int its_max, vec& eigenvals_A) ;

//////////////////  MAIN
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
        mat R(N, N, fill::eye);  // Matrix R to be filled with the eigenvectors
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
        A(3,4) = 51;
    std::cout << "\n Original Matrix A: \n";
    std::cout << A << "\n";

    /////////////
    // we want to solve Av = lambda v using arma::eig_sym
	vec eigval; // vector containing eigenvalues
	mat eigvec; // eigenvector
	eig_sym(eigval, eigvec, A);
    // eigenvalues in order
    std::cout << "Eigen-values: \n" << eigval << "\n";
    // eigevectors in a matrix
    std::cout << "Eigen-vectors: \n"<< eigvec << "\n";
    ////////////
    int p = 0;
    int q = 0;
    double max;
    max = offdiag(A, p, q,N);
    
    double tolerance = pow(10, -8); // pow(base, exponent)
    int its_max = 10000;
    vec eigenvals_A(N);//,sorted_eigenvals(N); // defining vector with N elements
    R = eigen_solver_loop(tolerance, A, R, N, its_max, eigenvals_A);

    

    // sorting_bysize(R,eigenvals_A, sorted_R, sorted_eigenvals_A);
    std::cout << "New matrix A: \n" << A << "\n";
    std::cout << "Matrix R: \n" << R << "\n";
    

////////////////////////////////
    // using armadillo to note the indices 
    uvec indices = sort_index(eigenvals_A);
    std::cout << "Sorted (now just index) eigenvals A: \n"<< indices << "\n";
    vec sorted_eigenvals(N);
    mat sorted_R(N,N);
    for (int i = 0; i<N; i++){
        sorted_eigenvals(i) = eigenvals_A(indices(i));
        for(int j=0; j<N;j++){
            sorted_R(i,j) = R(indices(i),j);
        }
    }
    std::cout << "Sorted matrix R \n" << sorted_R << "\n";
    ////////////////////////////////

    //A.print(std::cout);
    //std::cout << "jacobi_rotate matrix R (????idk) \n";
    //R.print(std::cout);
    return 0;
}


////////////////// OFF DIAGONAL
        // Determine the the max off-diagonal element of a symmetric matrix A
        // - Saves the matrix element indicies to k and l 
        // - Returns absolute value of A(k,l) as the function return value
    double offdiag(mat& A, int& p, int& q, int N) ///////////// OPPGAVE 3!!!   // & er reference/pointer
    {
            double max;
            for (int i = 0; i < N; ++i) // Looking at position A(0,N-1) to A(N-1,0)
            {
                for ( int j = i+1; j < N; ++j)
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

////////////////// JACOBI ROTATION
    // Performs a single Jacobi rotation, to "rotate away"
    // the off-diagonal element at A(k,l).
    // - Assumes symmetric matrix, so we only consider k < l
    // - Modifies the input matrices A and R
void jacobi_rotate(mat& A, mat& R, int k, int l, int N)
{   /*/ MATRIX:
        // a_kk         a_lk   (l > k)

        // a_kl         a_ll

             // når maxoffdiag < 10^-8 ish skal det stoppe, så printer vi matrisen
             // Telle antall ganger iterasjonene går for å få den roterte matrisen*/
                double s, c;
                if ( A(k,l) != 0.0 ) 
                {    
                    double t, tau; // t for tangent
                    tau = (A(l,l) - A(k,k))/(2*A(k,l));

                    if ( tau >= 0 )         
                        {t = 1.0/(tau + sqrt(1.0 + tau*tau));}   // se slutt av fl notat
                    else 
                        {t = -1.0/(-tau + sqrt(1.0 + tau*tau));}
                    
                    c = 1/sqrt(1+t*t);  // c for cosine
                    s = c*t;            // s for sine
                } 
                else 
                {
                    c = 1.0;    // c for cosine when A(k,l) = 0
                    s = 0.0;    // s for sine when A(k,l) = 0
                }
                
                // Updating
                double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
                a_kk = A(k,k);
                a_ll = A(l,l);
                A(k,k) = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll;   // since Akl=Alk
                A(l,l) = s*s*a_kk + 2.0*c*s*A(k,l) + c*c*a_ll;
                A(k,l) = 0.0; // hard-coding non-diagonal elements by hand
                A(l,k) = 0.0; // same here
                
                for ( int i = 0; i < N; i++ ) // For the elements that are not kk, ll, kl, lk
                {
                    if ( i != k && i != l ) 
                    {
                        a_ik = A(i,k);
                        a_il = A(i,l);
                        A(i,k) = c*a_ik - s*a_il;
                        A(k,i) = A(i,k);
                        A(i,l) = c*a_il + s*a_ik;
                        A(l,i) = A(i,l);
                    }   
                }

                for (int i = 0; i<N; i++)
                {
                    r_ik = R(i,k);
                    r_il = R(i,l);
                    R(i,k) = c*r_ik - s*r_il;
                    R(i,l) = c*r_il + s*r_ik;
                }
        //std::cout << "jacobi_rotate A??  \n" << A;
}



////////////////// EIGEN SOLVER
mat eigen_solver_loop(double tolerance, mat& A, mat& R, int N, int its_max, vec& eigenvals_A) 
{
    int p,q,its; // its = iterations (for counting how many interations the code have to do)
    double max;
    its = 0;
    max = offdiag(A,p,q, N);
    while (tolerance < max && its < its_max) 
    {
        // One rotation for every max value that is bigger than the tolerance
        // we set.
        jacobi_rotate(A, R, p, q, N);
        max = offdiag(A, p, q, N);
        its = its + 1;      // Counting iterations
    }
    // Extracting the eigenvalues from the rotated A
    //for (int i = 0; i < N; i++)
    //    {eigenvals_A(i) = A(i,i);}

    std::cout << "\n Iterations: " << its << "\n" << "\n";
    return R;
}