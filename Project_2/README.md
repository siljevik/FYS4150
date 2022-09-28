# Project 2 - FYS4150/3150

## Introduction
The main topic of this project is the scaling of equations, eigenvalue problems and code testing. Where we want to solve matrices for their eigenvalues, eigenvectors, and to later implement
our results in the Jacobi rotation algorithm (with code!).

All the problems (except problem 1) is saved, respectivly, in ../2, ../3, and so on. And in order to run the code, one has to compile the *.cpp file in the terminal or in a C++ compiler.

Coding-language(s): C++

## Problem 1
Problem one was solved analytically, and the solution is given in a seperate report handed in on canvas.

## Problem 2

In the ../2 directory; the implemented tridiagonal matrix and its eigenvalues/eigenvectors is saved as the main.cpp file. 

When running the code, the eigenvalues and -vectors will be printed out with the solutions' matrix. In this case; A.

The matrix- and eigenvalue-problems are solved using armadillo. Especially with the function; arma::eig_sym 

## Problem 3

In this problem we look at the off-diagonal of our matrix, and try to find the largest absolute value of the off-diagonal elements.

To check the code, we test it on a custom matrix B. (Given in the project description/taskpaper and defined as A)

## Problem 4

We implement Jacobi's rotational algorithm to solve the eigenvalue-problem (6), and check an 6x6 matrix A. The results should agree with the analytical results for N=6.

## Problem 5

Here we want to run the program for different matrices of size NxN, and from this we want to estimate how the number of required transformations scale with the size N.
We run the program main.cpp where we will get the required iterations needed to reach values below our set limit (pow(10,-8)). We select points for N= 5, 10, ..., 70, 75 and plot the 
iterations against the matrix size, N, in Python.


## Problem 6

We want to use the Jacobi code to solve eq. (6) for discretized x with n =10. The three eigenvectors corresponding to the three lowest eigenvalues are plotted with the vectorelements
v_i against the corresponding positions x_i. 

We include the boundrarypoints (x0,v0) and (xn,vn), and plot the corresponding analythical eigenvectors in the same plot.   (Not completed!)
