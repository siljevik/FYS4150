#############################
#                           #
#    POINT 2 IN PROBLEM 8   #
#                           #
#############################
import matplotlib.pyplot as plt

# The filenames of the files the data will be extracted from
filename_2 = "two_particle_with.txt" # With interactions
filename_1 = "two_particle_without.txt" # Without interactions



# This plots two particles movement in a Penning trap, in the 
# z-direction, as a function of time.
def two_particle_penning_xy(filename_1,filename_2):

    # Creating an empty lists for the x and y positions for both the
    # particles in the Penning trap with and without interactions
    x_w = []
    y_w = []
    x_wo = []
    y_wo = []
    


    ########################################################
    #                   WITH INTERACTIONS                  #
    # Extracts the data for the movement with interactions #
    ########################################################
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


    ########################################################
    #                 WITHOUT INTERACTIONS                 #
    # Extracts the data for the movement without           #
    # interactions                                         #
    ########################################################
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



    ######################################################
    #                   PLOTTING                         #
    # Plotting all 4 plots in a 2x2 grid in one svg file #
    ######################################################
    
    # 2 particles with interactions
    plt.plot(x_w, y_w, color='mediumslateblue')
    plt.title("Two-particle movement in the xy-plane, with interactions")
    plt.xlabel("x [\u03bcm]")
    plt.ylabel("y [\u03bcm]")
    plt.savefig("two_particle_xy_w.svg") # Saves it as a .svg-file
    plt.show()

    # 2 particles without interaction
    plt.plot(x_wo, y_wo, color='lightseagreen')
    plt.title("Two-particle movement in the xy-plane, without interactions")
    plt.xlabel("x [\u03bcm]")
    plt.ylabel("y [\u03bcm]")
    plt.savefig("two_particle_xy_wo.svg") # Saves it as a .svg-file
    plt.show()
    
    return x_w, y_w, x_wo, y_wo