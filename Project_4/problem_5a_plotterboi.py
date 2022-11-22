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
for line_1 in open(filename_1, 'r'):
    cols_1 = [i for i in line_1.split()]
    x_1.append(cols_1[0])
    y_1.append(int(cols_1[1]))
    z_1.append(int(cols_1[2]))

# For the initially ordered lattice for temperature T = 2.4 J/kB
for line_2 in open(filename_2, 'r'):
    cols_2 = [i for i in line_2.split()]
    x_2.append(cols_2[0])
    y_2.append(int(cols_2[1]))
    z_2.append(int(cols_2[2]))

# For the initially unordered lattice for temperature T = 1.0 J/kB
for line_3 in open(filename_3, 'r'):
    cols_3 = [i for i in line_3.split()]
    x_3.append(cols_3[0])
    y_3.append(int(cols_3[1]))
    z_3.append(int(cols_3[2]))

# For the initially unordered lattice for temperature T = 2.4 J/kB
for line_4 in open(filename_4, 'r'):
    cols_4 = [i for i in line_4.split()]
    x_4.append(cols_4[0])
    y_4.append(int(cols_4[1]))
    z_4.append(int(cols_4[2]))



########################################
#               Plotting               #
########################################

# PLOT 1: Energy development of an ordered lattice
plt.plot(np.asarray(x_1,float), np.asarray(y_1,float), color='orange', label='T = 1.0 J/kB')
plt.plot(np.asarray(x_2,float), np.asarray(y_2,float), color='sandybcoln', label='T = 2.4 J/kB')
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
plt.savefig("unrdered_epsilon.svg", format="svg")
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