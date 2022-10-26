#############################
#                           #
#    POINT 3 IN PROBLEM 8   #
#                           #
#############################
from platform import platform
import matplotlib.pyplot as plt
import numpy as np

# The filenames of the files the data will be extracted from
filename_1 = "two_particle_data_w_int.txt" # With interactions
filename_2 = "two_particle_data_no_int.txt" # Without interactions
# Both files has this header (the comment under is just to have 
# somewhat control of placement numbering):
# Time, x1, x2, y1, y2, z1, z2, vx1, vx2, vy1, vy2, vz1, vz2
#    0,  1,  2,  3,  4,  5,  6,   7,   8,   9,  10,  11,  12



# This plots two particles movement in a Penning trap, in the x,vx-plane and the z,vz-plane
def two_particle_penning_xzv(filename_1,filename_2):

    # Creating an empty lists for both the particles in the Penning trap with (w) and 
    # without (wo) interactions.
    # Particle 1:
    x1_w = []
    z1_w = []
    vx1_w = []
    vz1_w = []
    x1_wo = []
    z1_wo = []
    vx1_wo = []
    vz1_wo = []
    # Particle 2:
    x2_w = []
    z2_w = []
    vx2_w = []
    vz2_w = []
    x2_wo = []
    z2_wo = []
    vx2_wo = []
    vz2_wo = []
    


    ######################################################
    #                   WITH INTERACTIONS                #
    # Extracts the data for the trajectory with          #
    # interactions                                       #
    ######################################################
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
            z1_w_val = line_cont_ls[5]
            vx1_w_val = line_cont_ls[7]
            vz1_w_val = line_cont_ls[11]
            x2_w_val = line_cont_ls[2]
            z2_w_val = line_cont_ls[6]
            vx2_w_val = line_cont_ls[8]
            vz2_w_val = line_cont_ls[12]

            # Adding the values to the belonging lists
            x1_w.append(x1_w_val)
            z1_w.append(z1_w_val)
            vx1_w.append(vx1_w_val)
            vz1_w.append(vz1_w_val)
            x2_w.append(x2_w_val)
            z2_w.append(z2_w_val)
            vx2_w.append(vx2_w_val)
            vz2_w.append(vz2_w_val)



    ######################################################
    #                 WITHOUT INTERACTIONS               #
    # Extracts the data for the trajectory without       #
    # interactions                                       #
    ######################################################
    with open(filename_2) as f: 
        next(f)
        text = f.read()
        line_ls = text.split('\n')
        del line_ls[-1] # Last line is empty, so we delete it, bitch bye

        for j in line_ls:
            line_cont_ls = j.split(",") # List of values in the line

            # Extracting values from filename_2
            x1_wo_val = line_cont_ls[1]
            z1_wo_val = line_cont_ls[5]
            vx1_wo_val = line_cont_ls[7]
            vz1_wo_val = line_cont_ls[11]
            x2_wo_val = line_cont_ls[2]
            z2_wo_val = line_cont_ls[6]
            vx2_wo_val = line_cont_ls[8]
            vz2_wo_val = line_cont_ls[12]

            # Adding the values to the belonging lists
            x1_wo.append(x1_wo_val)
            z1_wo.append(z1_wo_val)
            vx1_wo.append(vx1_wo_val)
            vz1_wo.append(vz1_wo_val)
            x2_wo.append(x2_wo_val)
            z2_wo.append(z2_wo_val)
            vx2_wo.append(vx2_wo_val)
            vz2_wo.append(vz2_wo_val)

    # Returning all the complete lists
    return x1_w,z1_w,vx1_w,vz1_w,x1_wo,z1_wo,vx1_wo,vz1_wo,x2_w,z2_w,vx2_w,vz2_w,x2_wo,z2_wo,vx2_wo,vz2_wo





# RUNNING THE CODE
x1_w,z1_w,vx1_w,vz1_w,x1_wo,z1_wo,vx1_wo,vz1_wo,x2_w,z2_w,vx2_w,vz2_w,x2_wo,z2_wo,vx2_wo,vz2_wo = two_particle_penning_xzv(filename_1,filename_2)



######################################################
#                     PLOTTING                       #
# Plotting all 4 plots, saving them, one at a time   #
######################################################

# PLOT 1: For two particles' trajectories with interactions in the (x,vx) plane
plt.plot(np.asarray(x1_w,float), np.asarray(vx1_w,float), color='forestgreen', label='Particle 1')
plt.plot(np.asarray(x2_w,float), np.asarray(vx2_w,float), color='limegreen', label='Particle 2')
plt.title('Trajectories in the (x,$v_x$)-plane, with interactions')
plt.xlabel('x [\u03bcm]')
plt.ylabel('$v_x$ [\u03bcm/\u03bcs]')
plt.legend()
plt.savefig("two_particle_traj_x_w.svg", format="svg") 
plt.close()

# PLOT 2: For two particles' trajectories without interactions in the (x,vx) plane
plt.plot(np.asarray(x1_wo,float), np.asarray(vx1_wo,float), color='teal', label='Particle 1')
plt.plot(np.asarray(x2_wo,float), np.asarray(vx2_wo,float), color='mediumturquoise', label='Particle 2')
plt.title('Trajectories in the (x,$v_x$)-plane, without interactions')
plt.xlabel('x [\u03bcm]')
plt.ylabel('$v_x$ [\u03bcm/\u03bcs]')
plt.legend()
plt.savefig("two_particle_traj_x_wo.svg", format="svg") 
plt.close()

# PLOT 3: For two particles' trajectories with interactions in the (z,vz) plane
plt.plot(np.asarray(z1_w,float), np.asarray(vz1_w,float), color='salmon', label='Particle 1')
plt.plot(np.asarray(z2_w,float), np.asarray(vz2_w,float), color='red', label='Particle 2')
plt.title('Trajectories in the (z,$v_z$)-plane, with interactions')
plt.xlabel('z [\u03bcm]')
plt.ylabel('v\_{z}[\u03bcm/\u03bcs]')
plt.legend()
plt.savefig("two_particle_traj_z_w.svg", format="svg") 
plt.close()
     
# PLOT 4: For two particles' trajectories without interactions in the (z,vz) plane
plt.plot(np.asarray(z1_wo,float), np.asarray(vz1_wo,float), color='darkmagenta',label='Particle 1')
plt.plot(np.asarray(z2_wo,float), np.asarray(vz2_wo,float), color='orchid', label='Particle 2')
plt.title('Trajectories in the (x,$v_x$)-plane, without interactions')
plt.xlabel('z [\u03bcm]')
plt.ylabel('$v_z$[\u03bcm/\u03bcs]')
plt.legend()
plt.savefig("two_particle_traj_z_wo.svg", format="svg") 
plt.close()