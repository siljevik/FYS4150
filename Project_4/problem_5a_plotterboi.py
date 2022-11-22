######################################################
##          PROJECT 4 - PROBLEM 5a PLOTTER          ##
######################################################

## This Python code uses the .txt-files:
##   - equilibrium_time_T_1_0.txt
##   - equilibrium_time_T_2_4.txt
## or the txt-files
##   - random_time_T_1_0.txt
##   - random_time_T_1_4.txt
## where all files has these placements:
## col[0] = MC-cycles       col[1] =  <ϵ>     col[2] = <|m|>¨
# To create two plots: one with the numerical estimate <ϵ> with the number of MC-cycles (MC = Monte Carlo)

import matplotlib.pyplot as plt
import numpy as np
from turtle import color

# Creating empty lists to fill them up with data from the files

# 1 for T = 1.0 J/k_B
x_1 = [] # MC-cycles
y_1 = [] # <ϵ>
z_1 = [] # <|m|>

# 2 for T = 2.4 J/k:B
x_2 = [] # MC-cycles
y_2 = [] # <ϵ>
z_2 = [] # <|m|>

# 1 for T = 1.0 J/k_B
x_3 = [] # MC-cycles
y_3 = [] # <ϵ>
z_3 = [] # <|m|>

# 2 for T = 2.4 J/k:B
x_4 = [] # MC-cycles
y_4 = [] # <ϵ>
z_4 = [] # <|m|>



#################################
#           Data-files          #
#################################

# For lattice where all elements are initially ordered (has the same spin-direction)
filename_1 = 'equilibrium_time_T_1_0.txt'
filename_2 = 'equilibrium_time_T_2_4.txt'
# For lattice where all elements are intiially unordered (random spin-directions)
filename_3 = 'random_time_T_1_0.txt'
filename_4 = 'random_time_T_2_4.txt'



#################################################
#       Extracting information from files       #
#################################################

# For the initially ordered lattice for temperature T = 1.0 J/kB
with open(filename_1) as f_1: 
    next(f_1) # Skips the first line (since it is just a header)
    text_1 = f_1.read()
    # Splitting it by \n, to get a vector containing each line/row of txt file
    line_ls_1 = text_1.split('\n')
    del line_ls_1[-1] # The last line is empty, so we delete it
    for i_1 in line_ls_1:
        line_cont_ls_1 = i_1.split(" ") # List of values in the line
        x_1.append(line_cont_ls_1[0])
        y_1.append(line_cont_ls_1[1])
        z_1.append(line_cont_ls_1[2])

# For the initially ordered lattice for temperature T = 2.4 J/kB
with open(filename_2) as f_2: 
    next(f_2) 
    text_2 = f_2.read()
    line_ls_2 = text_2.split('\n')
    del line_ls_2[-1] 
    for i_2 in line_ls_2:
        line_cont_ls_2 = i_2.split(" ")
        x_2.append(line_cont_ls_2[0])
        y_2.append(line_cont_ls_2[1])
        z_2.append(line_cont_ls_2[2])

# For the initially unordered lattice for temperature T = 1.0 J/kB
with open(filename_3) as f_3: 
    next(f_3)
    text_3 = f_3.read()
    line_ls_3 = text_3.split('\n')
    del line_ls_3[-1]
    for i_3 in line_ls_3:
        line_cont_ls_3 = i_3.split(" ")
        x_3.append(line_cont_ls_3[0])
        y_3.append(line_cont_ls_3[1])
        z_3.append(line_cont_ls_3[2])

# For the initially unordered lattice for temperature T = 2.4 J/kB
with open(filename_4) as f_4: 
    next(f_4)
    text_4 = f_4.read()
    line_ls_4 = text_4.split('\n')
    del line_ls_4[-1]
    for i_4 in line_ls_4:
        line_cont_ls_4 = i_4.split(" ")
        x_4.append(line_cont_ls_4[0])
        y_4.append(line_cont_ls_4[1])
        z_4.append(line_cont_ls_4[2])



########################################
#               Plotting               #
########################################

# PLOT 1: Energy development of an ordered lattice
plt.plot(np.asarray(x_1,float), np.asarray(y_1,float), color='IndianRed', label='T = 1.0 J/kB')
plt.plot(np.asarray(x_2,float), np.asarray(y_2,float), color='LightSalmon', label='T = 2.4 J/kB')
plt.title('Energy development of an unordered lattice \n over Monte Carlo cycles')
plt.xlabel('Number of Monte Carlo cycles')
plt.ylabel('<ϵ>')
plt.legend()
plt.savefig("ordered_epsilon.svg", format="svg")
plt.close()

# PLOT 2: Magnetization development of an ordered lattice
plt.plot(np.asarray(x_1,float), np.asarray(z_1,float), color='palevioletred', label='T = 1.0 J/kB')
plt.plot(np.asarray(x_2,float), np.asarray(z_2,float), color='hotpink', label='T = 2.4 J/kB')
plt.title('Magnetization development of an ordered lattice \n over Monte Carlo cycles')
plt.xlabel('Number of Monte Carlo cycles')
plt.ylabel('<|m|>')
plt.legend()
plt.savefig("ordered_m.svg", format="svg") 
plt.close()

# PLOT 3: Energy development of an unordered lattice
plt.plot(np.asarray(x_3,float), np.asarray(y_3,float), color='forestgreen', label='T = 1.0 J/kB')
plt.plot(np.asarray(x_4,float), np.asarray(y_4,float), color='limegreen', label='T = 2.4 J/kB')
plt.title('Energy development of an unordered lattice \n over Monte Carlo cycles')
plt.xlabel('Number of Monte Carlo cycles')
plt.ylabel('<ϵ>')
plt.legend()
plt.savefig("unordered_epsilon.svg", format="svg")
plt.close()

# PLOT 4: Magnetization development of an unordered lattice
plt.plot(np.asarray(x_3,float), np.asarray(z_3,float), color='teal', label='T = 1.0 J/kB')
plt.plot(np.asarray(x_4,float), np.asarray(z_3,float), color='mediumturquoise', label='T = 2.4 J/kB')
plt.title('Magnetization development of an unordered lattice \n over Monte Carlo cycles')
plt.xlabel('Number of Monte Carlo cycles')
plt.ylabel('<|m|>')
plt.legend()
plt.savefig("unordered_m.svg", format="svg") 
plt.close()