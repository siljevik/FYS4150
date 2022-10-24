# POINT 1 IN PROBLEM 8
import matplotlib.pyplot as plt

# The filename of the file the data will be extracted from
filename = "single_particle_z.txt"

# This plots a single particles movement in a Penning trap, in the 
# z-direction, as a function of time.
def single_particle_penning_t_z(filename):

    # Creating an empty lists for the z-value and the time
    t = []
    z = []
    

    # Opens the file and reads it
    with open(filename, encoding="utf8") as f: 
        text = f.read()
        # Splitting it by \n
        line_ls = text.split('\n')

        for i in line_ls:
            line_cont_ls = i.split() # List of values in the line

            # Finding the t and z value from this line (position 1????)
            t_val = line_cont_ls[0]
            z_val = line_cont_ls[3]

            # Adding the t and z value to the list
            t.append(t_val)
            z.append(z_val)

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