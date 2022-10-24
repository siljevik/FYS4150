import matplotlib.pyplot as plt
import numpy as np


# This plots a single particles movement in a Penning trap, in the 
# z-direction, as a function of time.
def single_particle_penning_t_z(filename):

# Creating an empty lists for the z-value an the time
z = []

# for the time:
num_t = 0
tot_t = 50 # us

# Opens the file and reads it
with open(filename, encoding="utf8") as f: 
    text = f.read()
    # Splitting it by \n
    line_ls = text.split('\n')

    for i in line_ls:
        line_cont_ls = i.split() # List of values in the line

        # Finding the z value from this line (position 1????)
        z_val = line_cont_ls[1] #?????
        # Adding the z value to the list
        z.append(z_val)

        # Adding the to the time
        num_t += 1
    
    # Creating the list of times
    t = np.linspace(0,tot_t,num=num_t)


    # Plots the lines as solid green and dotted magneta
    plt.plot(t, z)
    # naming the plot and axes
    plt.title("Single particle movement")
    plt.xlabel("Time [\u03bcs]")
    plt.ylabel("Movement in z-direction[\u03bcm]")
    # Saves and shows the plot
    plt.savefig("single_particle_zt.png") 
    plt.show()

return z