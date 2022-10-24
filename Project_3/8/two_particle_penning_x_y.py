# POINT 2 IN PROBLEM 8
import matplotlib.pyplot as plt

# The filenames of the files the data will be extracted from
filename_2 = "two_particle_with.txt" # With interactions
filename_1 = "two_particle_without.txt" # Without interactions


# This plots two particles movement in a Penning trap, in the 
# z-direction, as a function of time.
def single_particle_penning_t_z(filename_1,filename_2):

    # Creating an empty lists for the x and y positions for both the
    # particles in the Penning trap with and without interactions
    x_w = []
    y_w = []
    x_wo = []
    y_wo = []
    

    # Opens the filename_1 and reads it
    with open(filename_1, encoding="utf8") as f: 
        text = f.read()
        # Splitting it by \n
        line_ls = text.split('\n')

        for i in line_ls:
            line_cont_ls = i.split() # List of values in the line

            # Finding the t and z value from this line (position 1????)
            x_w_val = line_cont_ls[0]
            y_w_val = line_cont_ls[3]

            # Adding the t and z value to the list
            x_w.append(x_w_val)
            y_w.append(y_w_val)

        # Plots the lines as solid green and dotted magneta
        plt.plot(x_w, y_w)
        # naming the plot and axes
        plt.title("Two-particle movement in the xy-plane, with interactions")
        plt.xlabel("x [\u03bcm]")
        plt.ylabel("y [\u03bcm]")
        # Saves and shows the plot
        plt.savefig("two_particle_xy_w.svg") 
        plt.show()


    # Opens the file and reads it
    with open(filename_2, encoding="utf8") as f: 
        text = f.read()
        # Splitting it by \n
        line_ls = text.split('\n')

        for i in line_ls:
            line_cont_ls = i.split() # List of values in the line

            # Finding the t and z value from this line (position 1????)
            x_wo_val = line_cont_ls[0]
            y_wo_val = line_cont_ls[3]

            # Adding the t and z value to the list
            x_wo.append(x_wo_val)
            y_wo.append(y_wo_val)

        # Plots the lines as solid green and dotted magneta
        plt.plot(x_wo, y_wo)
        # naming the plot and axes
        plt.title("Two-particle movement in the xy-plane, without interactions")
        plt.xlabel("x [\u03bcm]")
        plt.ylabel("y [\u03bcm]")
        # Saves and shows the plot
        plt.savefig("two_particle_xy_wo.svg") 
        plt.show()

    return x_w, y_w, x_wo, y_wo