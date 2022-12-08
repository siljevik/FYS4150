#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <random>
#include <vector>
#include <map>


#include "header.hpp"


arma::vec Header::indextranslator(int M){

        // defining an empty vector
        std::vector<double> un_vec;

        // fill the vector with vectors
        for(int k = 0; k < M-2; k++)
        {

		int uij = 1;
                for(int i = 0; i < M; i++)
                {       for(int j = 0; j < M; j++)
                        {
                                //append uij values into the uij_vec vector which will be appended into un_vec
                                un_vec.push_back(uij);
                        }//end of j-loop
                }// end of i-loop
        }// end of k-loop

        return un_vec;
}
