#############################
#                           #
#    POINT 2 IN PROBLEM 8   #
#                           #
#############################
import matplotlib.pyplot as plt
import numpy as np

# The filenames of the files the data will be extracted from
filename_1 = "two_particle_data_w_int.txt" # With interactions
filename_2 = "two_particle_data_no_int.txt" # Without interactions
# Both files has this header (the comment under is just to have 
# somewhat control of placement numbering):
# Time, x1, x2, y1, y2, z1, z2, vx1, vx2, vy1, vy2, vz1, vz2
#    0,  1,  2,  3,  4,  5,  6,   7,   8,   9,  10,  11,  12



# This plots two particles movement in a Penning trap in the xy-plane
def two_particle_penning_xy(filename_1,filename_2):

    # Creating an empty lists for the x and y positions for both the
    # particles in the Penning trap with and without interactions
    # Particle 1:
    x1_w = []
    y1_w = []
    x1_wo = []
    y1_wo = []
    # Particle 2:
    x2_w = []
    y2_w = []
    x2_wo = []
    y2_wo = []


    ########################################################
    #                   WITH INTERACTIONS                  #
    # Extracts the data for the movement with interactions #
    ########################################################
    with open(filename_1) as f: 
        next(f) # Skipping the header (top) row
        text = f.read()

        # Splitting it by \n
        line_ls = text.split('\n')
        del line_ls[-1] # Last line is empty, so we delete it, bitch bye

        for i in line_ls:
            line_cont_ls = i.split(",") # List of values in the line

            # Extracting values from filename_1
            x1_w_val = line_cont_ls[1]
            y1_w_val = line_cont_ls[3]
            x2_w_val = line_cont_ls[2]
            y2_w_val = line_cont_ls[4]

            # Adding the values to the belonging lists
            x1_w.append(x1_w_val)
            y1_w.append(y1_w_val)
            x2_w.append(x2_w_val)
            y2_w.append(y2_w_val)


    ########################################################
    #                 WITHOUT INTERACTIONS                 #
    # Extracts the data for the movement without           #
    # interactions                                         #
    ########################################################
    # Opens the file and reads it
    with open(filename_2) as f:
        next(f) 
        text = f.read()
        # Splitting it by \n
        line_ls = text.split('\n')
        del line_ls[-1] # Last line is empty, so we delete it, bitch bye

        for j in line_ls:
            line_cont_ls = j.split(",") # List of values in the line

            # Finding the t and z value from this line (position 1????)
            x1_wo_val = line_cont_ls[1]
            y1_wo_val = line_cont_ls[3]
            x2_wo_val = line_cont_ls[2]
            y2_wo_val = line_cont_ls[4]

            # Adding the t and z value to the list
            x1_wo.append(x1_wo_val)
            y1_wo.append(y1_wo_val)
            x2_wo.append(x2_wo_val)
            y2_wo.append(y2_wo_val)

    
    return x1_w, y1_w, x2_w, y2_w, x1_wo, y1_wo, x2_wo, y2_wo


# RUNNING THE CODE
x1_w, y1_w, x2_w, y2_w, x1_wo, y1_wo, x2_wo, y2_wo = two_particle_penning_xy(filename_1,filename_2)

######################################################
#                     PLOTTING                       #
# Plotting both plots, saving them, one at a time    #
######################################################

# PLOT 1: For two particles' trajectories with interactions in the (x,vx) plane
plt.plot(np.asarray(x1_w,float), np.asarray(y1_w,float), color='orange', label='Particle 1')
plt.plot(np.asarray(x2_w,float), np.asarray(y2_w,float), color='sandybrown', label='Particle 2')
plt.title('Trajectories in the xy-plane, with interactions')
plt.xlabel('x [\u03bcm]')
plt.ylabel('y [\u03bcm]')
plt.legend()
plt.savefig("two_particle_traj_xy_w.svg", format="svg") 
plt.close()

# PLOT 2: For two particles' trajectories without interactions in the (x,vx) plane
plt.plot(np.asarray(x1_wo,float), np.asarray(y1_wo,float), color='palevioletred', label='Particle 1')
plt.plot(np.asarray(x2_wo,float), np.asarray(y2_wo,float), color='hotpink', label='Particle 2')
plt.title('Trajectories in the xy-plane, without interactions')
plt.xlabel('x [\u03bcm]')
plt.ylabel('y [\u03bcm]')
plt.legend()
plt.savefig("two_particle_traj_xy_wo.svg", format="svg") 
plt.close()