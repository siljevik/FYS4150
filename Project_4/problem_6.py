#############################################
##          PROJECT 4 - PROBLEM 6          ##
#############################################
import matplotlib.pyplot as plt
import numpy as np
from turtle import color

# Creating empty lists to fill them up with data from the files
# 1 for T = 1.0 J/k_B
y_1 = [] # <ϵ>
# 2 for T = 2.4 J/k:B
y_2 = [] # <ϵ>
# 1 for T = 1.0 J/k_B
y_3 = [] # <ϵ>
# 2 for T = 2.4 J/k:B
y_4 = [] # <ϵ>



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
        y_1.append(line_cont_ls_1[1])

# For the initially ordered lattice for temperature T = 2.4 J/kB
with open(filename_2) as f_2: 
    next(f_2) 
    text_2 = f_2.read()
    line_ls_2 = text_2.split('\n')
    del line_ls_2[-1] 
    for i_2 in line_ls_2:
        line_cont_ls_2 = i_2.split(" ")
        y_2.append(line_cont_ls_2[1])

# For the initially unordered lattice for temperature T = 1.0 J/kB
with open(filename_3) as f_3: 
    next(f_3)
    text_3 = f_3.read()
    line_ls_3 = text_3.split('\n')
    del line_ls_3[-1]
    for i_3 in line_ls_3:
        line_cont_ls_3 = i_3.split(" ")
        y_3.append(line_cont_ls_3[1])

# For the initially unordered lattice for temperature T = 2.4 J/kB
with open(filename_4) as f_4: 
    next(f_4)
    text_4 = f_4.read()
    line_ls_4 = text_4.split('\n')
    del line_ls_4[-1]
    for i_4 in line_ls_4:
        line_cont_ls_4 = i_4.split(" ")
        y_4.append(line_cont_ls_4[1])



########################################
#               Plotting               #
########################################

# PLOT 1: Energy histogram with T=1.0, for initially unordered lattice
plt.hist(np.asarray(y_1,float),bins=100,color='IndianRed', label='T = 1.0 J/kB')
plt.title('probability distribution of ϵ (for intially ordered lattice)')
plt.xlabel('<ϵ> [J]')
plt.ylabel('# of ϵ')
plt.legend()
plt.savefig("prob6_T_1_ord.svg", format="svg")
plt.close()

# PLOT 2: Energy histogram with T=2.4, for initially ordered lattice
plt.hist(np.asarray(y_2,float),bins=100,color='hotpink', label='T = 2.4 J/kB')
plt.title('probability distribution of ϵ (for intially ordered lattice)')
plt.xlabel('<ϵ> [J]')
plt.ylabel('# of ϵ')
plt.legend()
plt.savefig("prob6_T_2_4_ord.svg", format="svg")
plt.close()

# PLOT 3: Energy development of an unordered lattice
plt.hist(np.asarray(y_3,float),bins=100,color='forestgreen', label='T = 1.0 J/kB')
plt.title('probability distribution of ϵ (for intially unordered lattice)')
plt.xlabel('<ϵ> [J]')
plt.ylabel('# of ϵ')
plt.legend()
plt.savefig("prob6_T_1_unord.svg", format="svg")
plt.close()

# PLOT 4: Magnetization development of an unordered lattice
plt.hist(np.asarray(y_4,float),bins=100,color='teal', label='T = 2.4 J/kB')
plt.title('probability distribution of ϵ (for intially unordered lattice)')
plt.xlabel('<ϵ> [J]')
plt.ylabel('# of ϵ')
plt.legend()
plt.savefig("prob6_T_2_4_unord.svg", format="svg")
plt.close()