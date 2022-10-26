#############################
#                           #
#    POINT 1 IN PROBLEM 8   #
#                           #
#############################
from turtle import color
import matplotlib.pyplot as plt
import numpy as np

# The filename of the file the data will be extracted from
# made with single_test.cpp
filename = "single_particle_z.txt"
# The file header looks like this:
# Time, x_num, y_num, z_num, x_ana, y_ana, z_ana
# 0   , 1    , 2    , 3    , 4    , 5    , 6        <-- Column number

# This plots a single particles movement in a Penning trap, in the 
# z-direction, as a function of time.
def single_particle_penning_t_z(filename):

    # Creating an empty lists for the z-value and the time
    t = []
    z = []
    zr = []

    
    

    # Opens the file and reads it
    with open(filename) as f: 
        next(f) # Skips the first line (since it is just a header)
        
        text = f.read()
        # Splitting it by \n, to get a vector containing each line/row of txt file
        line_ls = text.split('\n')
        del line_ls[-1] # Nothing to see here. Bitch bye

        for i in line_ls:
            line_cont_ls = i.split(",") # List of values in the line
            # Finding the t and z value from this line
            t_val = line_cont_ls[0]
            z_val = line_cont_ls[6]

            # Adding the t and z value to the list
            t.append(t_val)
            z.append(z_val)
    wz_sq = (2*2410000)/(40.078*(500**2))
    wz = np.sqrt(wz_sq)

    for j in t:
        z_val = 20*np.cos(wz * j)
        z.append(z_val)

    return t,z, zr

# For running the code:
t,z, zr = single_particle_penning_t_z(filename)

# Plots the lines
plt.plot(np.asarray(t, float), np.asarray(zr, float), color='orchid')
# naming the plot and axes
plt.title("Single particle movement")
plt.xlabel("Time [\u03bcs]")
plt.ylabel("Length in z-direction [\u03bcm]")

# Saves and shows the plot
plt.savefig("single_particle_zt.svg", format="svg") 
plt.show()