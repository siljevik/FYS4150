# Project 3 - Penning Trap

## Introduction


The purpouse of this project is to study the effects of the Penning trap and develop an object-oriented code to simulate a set of particles with charges(q) and masses(m) in the trap.

Task 1 to 4 is done analytically and answers are written down in the report, while 5 to 9 is the simulation.

### Problem 5
For the particle.hpp headerfile, the member variables q,m,vec r, and vec v are stored under the particle class, where the constructor assigns said values.

### Problem 6

PenningTrap class is created in penningtrap.hpp and contain the member variables B0, V0, d
and a vector to contain all the Particle objects in the Penning trap. Furthermore, the class
contains the memberfunctions for
        - External electric fiels
        - External magnetic field
        - The force due to the interaction among the particles

### Problem 7

Expanding the PennigTrap classfile with Runge-Kutta 4th order and forward Euler method to evolve the Penning trap simulation. These are defined as member functions in pennigtrap.hpp and directly implemented in the penningtrap.cpp


### Problem 8

Simulation of a single particle, and two particles with and without interactions. Both single_test.cpp and double_test.cpp includes files from header/ and func/

#### For the single particle

Compile w/: g++ single_test.cpp func/particle.cpp func/penningtrap.cpp func/analytical.cpp -std=c++11 -I include -o single_test.exe -larmadillo

The simulation of the single particle is conducted with the single_test.cpp file and executed w/ ./single_test.exe > data_single_text.txt.

#### for a pair of particles

Compile w/ : g++ double_test.cpp func/particle.cpp func/penningtrap.cpp func/analytical.cpp -std=c++11 -I include -o double_test.exe -larmadillo

To execute file and save data (without interactions): ./double_test.exe > two_particle_data_no_int.txt

To execute file and save data (with interactions): ./double_test.exe > two_particle_data_w_int.txt
#### Single particle w/ more steps

Compile w/: g++ r_err.cpp func/particle.cpp func/penningtrap.cpp func/analytical.cpp -std=c++11 -I include -o r_err.exe -larmadillo

execute w/: ./r_err.exe

By using the simulation results for the four different stepsizes, the error convergence is estimated for RK4 and FE, and printed out in the terminal.
(Should return the values FE = 1.443 and RK4 = 1)

* Lack graphs in 4 subplots to illustrate the relative error*


### Problem 9

The code for problem 9 is given in main.cpp  and is compiled and ran w/

- Compile: g++ -O2 main.cpp func/particle.cpp  -std=c++11 -I include -o main.exe -larmadillo
- Run: ./main.exe

We have the code ready and it runs/produces data, but no graphs to illustrate results...

## Folder 8

### single_particle_z.txt
Data file made with single_test.cpp (in the Proj_3_copy).

### two_particle_data_no_int.txt
Data file made with double_test.cpp (in the Proj_3_copy), for no interactions between particles.

### two_particle_data_w_int.txt
Data file made with double_test.cpp (in the Proj_3_copy), for with interactions between particles.

### single_particle_penning_t_z.py
Python code that plots for one particle in the penning trap with the single_particle_z.txt file.

### two_particle_penning_x_y.py
Plots the two particles' movement in the xy-plane, with and without interactions. Data extracted from two_particle_data_no_int.txt and two_particle_data_w_int.txt.

### two_particles_penning_traj.py
Plots the trajectories for two particles in the (x,v_x)- and (z,v_z)-planes, both with and without interactions. Data extracted from two_particle_data_no_int.txt and two_particle_data_w_int.txt.

### two_particle_3D.py
Python plot that plots the motion of two particles in the 3D-space, both with and wihtout interactions. Data extracted from two_particle_data_no_int.txt and two_particle_data_w_int.txt.

### All .svg-files
All of the .svg files are outputs from the Pyhton codes.
