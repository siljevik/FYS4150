# Project 4 - The Ising model in two Dimensions

The topic of this project is the Ising model in two dimensions. We will use this simple model to explore temperature-dependent behaviour in ferromagnets. A particular goal is to numerically estimate the critical temperature at which our 2D system undergoes a phase transition, from a magnetized phase to a phase with no net magnetization.

### Project focus

- The Markov Chain Monte Carlo method for sampling from probability distributions of many variables, and how the resulting samples can be used to approximate probability distributions and expectation values for derived quantities.
- Using parallelization to speed up our code.


## ~*~ C++ Code ~*~ 

### main.cpp
Compile with: g++ -O2 main.cpp -std=c++11 -I include -o main.exe -larmadillo

Run with: ./main.exe      

### MCMC_spin.hpp
Header file for the Markov Chain Monte Carlo simulations. 


### MCMC_spin.cpp
Function file for the Markov Chain Monte Carlo simulations.


### analytical.hpp
Header file for the analytical calculations.


### analytical.cpp
Function file for the analytical calculations.

### rand_state.cpp and main_omp_outer_loop.cpp
Attempt too parallelize the code, but didn't finish due to time constraints.


## ~*~ Plotting ~*~

### problem_5a_plotterboi.py
This Python code uses the .txt-files:
- equilibrium_time_T_1_0.txt
- equilibrium_time_T_2_4.txt
- random_time_T_1_0.txt
- random_time_T_1_4.txt
where all files has these placements:
col[0] = MC-cycles       col[1] =  <ϵ>     col[2] = <|m|>

To create two plots:
- The numerical estimate <ϵ> with the number of MC-cycles (MC = Monte Carlo)
- The magnetization <|m|> with the number of MC-cycles

### problem_6.py
Creates histograms based on the same .txt files as for problem_5a_plotterboi.py of the generated ϵ samples, using bin width 100.

### Data (.txt) files
All the .txt files are data-files created with main.cpp, with the content described under problem_5a_plotterboi.py above.
- random_time_T_2_4.txt
- random_time_T_1_0.txt
- equilibrium_time_T_2_4.txt
- equilibrium_time_T_1_0.txt

### Figure (.svg) files
the .svg files created in problem_5a_plotterboi:
- ordered_epsilon.svg
- ordered_m.svg
- unordered_epsilon.svg
- unordered_m.svg


The .svg files created in 
- prob6_T_1_ord.svg
- prob6_T_2_4_ord.svg
- prob6_T_1_unord.svg
- prob6_T_2_4_unord.svg

