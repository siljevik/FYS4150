///////////////////////////////////////////
//  FOR SILJES EYES ONLY:                //
//  Compile with:                        //
//  g++ -std=c++11 prob4.cpp -larmadillo //
//  Run with:                            //
//  ./a.out                              //
///////////////////////////////////////////

#include <iostream>
#include <armadillo>
using namespace arma;   // So we don't have to use 'arma::' before every Armadillo operation

// Declarations before our main
double offdiag(mat& A, int& p, int& q, int N); // The '&' denotes reference to an object, so we can change it
void jacobi_rotate(mat& A, mat& R, int k, int l, int N);
mat eigen_solver_loop(double tolerance, mat& A, mat& R, int N, int its_max, vec& eigenvals_A) ;



//////////////////////////////  MAIN   //////////////////////////////
int main() 
{
    int N = 6;          // Size of matrix A is NxN
    double x_min = 0.0; // min point
    double x_max = 1.0; // max point
    int  n = N + 1;
    double dx = x_max - x_min;
    double h = dx*1./n;
    mat A(N, N, fill::zeros); // Creating an empty matrix A, of size NxN, with Armadillo
    mat R(N, N, fill::eye);  // Matrix R to be filled with the eigenvectors
    // Creating two loops to fill the tridiagonal matrix with values:
    for (int i=0; i<N; i++)
        {A(i, i) =  2./(h*h);}
    // values,c ,diag.
    for (int i=0; i<(N-1) ; i++)
        {A(i+1, i) = -1./(h*h);
         A(i, i+1) = -1./(h*h);}
    // Just to have a look at the matrix we just made
    std::cout << "\n Original Matrix A: \n";
    std::cout << A << "\n";

    // Solving A*v = lambda*v using arma::eig_sym
	vec eigval; // vector containing eigenvalues
	mat eigvec; // matrix containint eigenvectors
	eig_sym(eigval, eigvec, A);
    // eigenvalues of our matrix, in increasing order
    std::cout << "Original Eigen-values: \n" << eigval << "\n";
    // Eigevectors of our matrix, where the first column belongs to first entry in the eigval
    // vector, and the second column belongs to the second entry etc.
    std::cout << "Original Eigen-vectors: \n"<< eigvec << "\n";
    
    int p = 0;
    int q = 0;
    double max;
    max = offdiag(A, p, q,N); // Finds the largest number/entry in the matrix that is not in the diagonal
    
    double tolerance = pow(10, -8); // pow(base, exponent) 'ten in the power of minus eight'
    int its_max = 10000;            // Let's not let our program go too crazy:D
    vec eigenvals_A(N);             // defining vector with N elements
    R = eigen_solver_loop(tolerance, A, R, N, its_max, eigenvals_A); // Creating a matrix where each column is an eigenvector

    

    // sorting_bysize(R,eigenvals_A, sorted_R, sorted_eigenvals_A);

    // Printing our updated matrixes neatly (imo)
    std::cout << "New matrix A: \n" << A << "\n";
    std::cout << "Matrix R: \n" << R << "\n";
    

    ////////////////////////////////
    uvec indices = sort_index(eigenvals_A); // using armadillo to note the indices 
    std::cout << "Sorted (now just index) eigenvals A: \n"<< indices << "\n";
    vec sorted_eigenvals(N);    // Creating a vector, sorted_eigenvals, of length N
    mat sorted_R(N,N);          // Creating a matrix, sorted_R, of size NxN
    for (int i = 0; i<N; i++){
        sorted_eigenvals(i) = eigenvals_A(indices(i));
        for(int j=0; j<N;j++){
            sorted_R(i,j) = R(indices(i),j);
        }
    }
    std::cout << "Sorted matrix R \n" << sorted_R << "\n";
    ////////////////////////////////
    return 0;
}





////////////////////////////// OFF DIAGONAL  //////////////////////////////
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
            double aij = fabs(A(i,j)); // Absolute value of A(i,j)
            if ( aij > max)
            {
                max = aij; p = i; q = j; // Updating our values/numbers
            }
        }
            
    }
    return max; // Returns the max off-diagonal element of A
}



//////////////////////////////  JACOBI ROTATION //////////////////////////////
// Performs a single Jacobi rotation, to "rotate away"
// the off-diagonal element at A(k,l).
// - Assumes symmetric matrix, so we only consider k < l
// - Modifies the input matrices A and R
void jacobi_rotate(mat& A, mat& R, int k, int l, int N)
{  
    double s, c;
    if ( A(k,l) != 0.0 ) 
    {    
        double t, tau; // t for tangent, tau for.. tau? haha
        tau = (A(l,l) - A(k,k))/(2*A(k,l));

        if ( tau >= 0 )         
            {t = 1.0/(tau + sqrt(1.0 + tau*tau));}   // from our calculations
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
    A(k,k) = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll;   // 2* since A(k,l)=A(l,k)
    A(l,l) = s*s*a_kk + 2.0*c*s*A(k,l) + c*c*a_ll;
    A(k,l) = 0.0; // Manually setting the non-diagonal values to 0
    A(l,k) = 0.0; 
                
    for ( int i = 0; i < N; i++ ) // For the elements that are not kk, ll, kl, lk
    {
        if ( i != k && i != l ) 
        {
            // Updating the non-kk/ll/kl/lk positions
            a_ik = A(i,k);
            a_il = A(i,l);
            A(i,k) = c*a_ik - s*a_il;
            A(k,i) = A(i,k);
            A(i,l) = c*a_il + s*a_ik;
            A(l,i) = A(i,l);
        }   
    }

    for (int i = 0; i < N; i++)
    {
        // Updating the values in the R matrix (each column will be an eigenvector of A)
        r_ik = R(i,k);
        r_il = R(i,l);
        R(i,k) = c*r_ik - s*r_il;
        R(i,l) = c*r_il + s*r_ik;
    }
}



////////////////////////////// EIGEN SOLVER //////////////////////////////
mat eigen_solver_loop(double tolerance, mat& A, mat& R, int N, int its_max, vec& eigenvals_A) 
{
    int p,q,its; // its = iterations (for counting how many interations the code have to do)
    double max;
    its = 0;
    max = offdiag(A,p,q, N);
    while (tolerance < max && its < its_max) 
    {
        // One rotation for every max value that is bigger than the tolerance we set.
        jacobi_rotate(A, R, p, q, N);
        max = offdiag(A, p, q, N);
        its = its + 1;      // Counting iterations
    }
    // Printing how many iterations our code had to do
    std::cout << "\n Iterations: " << its << "\n" << "\n";
    return R;
}