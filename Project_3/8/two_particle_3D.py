#############################
#                           #
#    POINT 4 IN PROBLEM 8   #
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



# This plots two particles movement in a Penning trap, in the 
# z-direction, as a function of time.
def two_particle_penning_3D(filename_1,filename_2):
    # Creating an empty lists for the x and y positions for both the
    # particles in the Penning trap with and without interactions
    # Particle 1:
    x1_w = []
    y1_w = []
    z1_w = []
    x1_wo = []
    y1_wo = []
    z1_wo = []
    # Particle 2:
    x2_w = []
    y2_w = []
    z2_w = []
    x2_wo = []
    y2_wo = []
    z2_wo = []



    ########################################################
    #                   WITH INTERACTIONS                  #
    # Extracts the data for the movement with interactions #
    ########################################################
    with open(filename_1) as f:
        next(f) 
        text = f.read()

        # Splitting it by \n
        line_ls = text.split('\n')
        del line_ls[-1] # Last line is empty, so we delete it, bitch bye

        for i in line_ls:
            line_cont_ls = i.split(",") # List of values in the line

            # Finding the t and z value from this line (position 1????)
            x1_w_val = line_cont_ls[1]
            y1_w_val = line_cont_ls[3]
            z1_w_val = line_cont_ls[5]
            x2_w_val = line_cont_ls[2]
            y2_w_val = line_cont_ls[4]
            z2_w_val = line_cont_ls[6]

            # Adding the t and z value to the list
            x1_w.append(x1_w_val)
            y1_w.append(y1_w_val)
            z1_w.append(z1_w_val)
            x2_w.append(x2_w_val)
            y2_w.append(y2_w_val)
            z2_w.append(z2_w_val)



    ########################################################
    #                 WITHOUT INTERACTIONS                 #
    # Extracts the data for the movement without           #
    # interactions                                         #
    ########################################################
    with open(filename_2) as f: 
        next(f) # Skipping the header (top) row
        text = f.read()

        # Splitting it by \n
        line_ls = text.split('\n')
        del line_ls[-1] # Last line is empty, so we delete it, bitch bye

        for j in line_ls:
            line_cont_ls = j.split(",") # List of values in the line


            # Extracting values from filename_2
            x1_wo_val = line_cont_ls[1]
            y1_wo_val = line_cont_ls[3]
            z1_wo_val = line_cont_ls[5]
            x2_wo_val = line_cont_ls[2]
            y2_wo_val = line_cont_ls[4]
            z2_wo_val = line_cont_ls[6]

            # Adding the values to the belonging lists
            x1_wo.append(x1_wo_val)
            y1_wo.append(y1_wo_val)
            z1_wo.append(z1_wo_val)
            x2_wo.append(x2_wo_val)
            y2_wo.append(y2_wo_val)
            z2_wo.append(z2_wo_val)

    
    return x1_w,y1_w,z1_w,x1_wo,y1_wo,z1_wo,x2_w,y2_w,z2_w,x2_wo,y2_wo,z2_wo



# RUNNING THE CODE
x1_w,y1_w,z1_w,x1_wo,y1_wo,z1_wo,x2_w,y2_w,z2_w,x2_wo,y2_wo,z2_wo = two_particle_penning_3D(filename_1,filename_2)

######################################################
#                   PLOTTING                         #
# Plotting all 4 plots in a 2x2 grid in one svg file #
######################################################

# PLOT 1: For two particles' trajectories with interactions in the 3D-plane
ax1 = plt.figure().add_subplot(projection='3d')
ax1.plot(np.asarray(x1_w,float), np.asarray(y1_w,float), np.asarray(z1_w,float), color='deepskyblue', label='Particle 1')
ax1.plot(np.asarray(x2_w,float), np.asarray(y2_w,float), np.asarray(z2_w,float), color='steelblue', label='Particle 2')
ax1.set_title('Trajectories in the 3D-plane, with interactions')
ax1.set_xlabel('x [\u03bcm]')
ax1.set_ylabel('y [\u03bcm]')
ax1.set_zlabel('z [\u03bcm]')
ax1.legend()
plt.savefig("two_particle_3D_w.svg", format="svg") 
plt.close()


# PLOT 1: For two particles' trajectories without interactions in the 3D-plane
ax2 = plt.figure().add_subplot(projection='3d')
ax2.plot(np.asarray(x1_wo,float), np.asarray(y1_wo,float), np.asarray(z1_wo,float), color='peru', label='Particle 1')
ax2.plot(np.asarray(x2_wo,float), np.asarray(y2_wo,float), np.asarray(z2_wo,float), color='saddlebrown', label='Particle 2')
ax2.set_title('Trajectories in the 3D-plane, without interactions')
ax2.set_xlabel('x [\u03bcm]')
ax2.set_ylabel('y [\u03bcm]')
ax2.set_zlabel('z [\u03bcm]')
ax2.legend()
plt.savefig("two_particle_3D_wo.svg", format="svg") 
plt.close()