#############################
#                           #
#    POINT 3 IN PROBLEM 8   #
#                           #
#############################
import matplotlib.pyplot as plt

# The filenames of the files the data will be extracted from
filename_2 = "two_particle_with_traj.txt" # With interactions
filename_1 = "two_particle_without_traj.txt" # Without interactions



# This plots two particles movement in a Penning trap, in the 
# z-direction, as a function of time.
def single_particle_penning_t_z(filename_1,filename_2):

    # Creating an empty lists for the x and y positions for both the
    # particles in the Penning trap with and without interactions
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
    with open(filename_1, encoding="utf8") as f: 
        text = f.read()
        # Splitting it by \n
        line_ls = text.split('\n')

        for i in line_ls:
            line_cont_ls = i.split() # List of values in the line

            # Finding the t and z value from this line (position 1????)
            x1_w_val = line_cont_ls[0]
            z1_w_val = line_cont_ls[3]
            vx1_w_val = line_cont_ls[0]
            vz1_w_val = line_cont_ls[3]
            x2_w_val = line_cont_ls[0]
            z2_w_val = line_cont_ls[3]
            vx2_w_val = line_cont_ls[0]
            vz2_w_val = line_cont_ls[3]

            # Adding the t and z value to the list
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
    with open(filename_2, encoding="utf8") as f: 
        text = f.read()
        # Splitting it by \n
        line_ls = text.split('\n')

        for i in line_ls:
            line_cont_ls = i.split() # List of values in the line

            # Finding the t and z value from this line (position 1????)
            x1_wo_val = line_cont_ls[0]
            z1_wo_val = line_cont_ls[3]
            vx1_wo_val = line_cont_ls[0]
            vz1_wo_val = line_cont_ls[3]
            x2_wo_val = line_cont_ls[0]
            z2_wo_val = line_cont_ls[3]
            vx2_wo_val = line_cont_ls[0]
            vz2_wo_val = line_cont_ls[3]

            # Adding the t and z value to the list
            x1_wo.append(x1_wo_val)
            z1_wo.append(z1_wo_val)
            vx1_wo.append(vx1_wo_val)
            vz1_wo.append(vz1_wo_val)
            x2_wo.append(x2_wo_val)
            z2_wo.append(z2_wo_val)
            vx2_wo.append(vx2_wo_val)
            vz2_wo.append(vz2_wo_val)



    ######################################################
    #                   PLOTTING                         #
    # Plotting all 4 plots in a 2x2 grid in one svg file #
    ######################################################
    
    fig, axs = plt.subplots(2, 2)
    # Giving the whole cllection of plots a title
    fig.suptitle('Trajectories of the two particles...')

    # Plot for two particles' trajectories with interactions in the (x,vx) plane
    axs[0, 0].plot(x1_w, vx1_w, color='forestgreen', label='Particle 1')
    axs[0, 0].plot(x2_w, vx2_w, color='limegreen', label='Particle 2')
    axs[0, 0].set_title('... in the (x,$v\_{x}$) plane, with interactions')
    axs[0, 0].xlabel('x [\u03bcm]')
    axs[0, 0].ylabel('$v\_{x}$ [\u03bcm/\u03bcms]')
    axs[0, 0].legend()

    # Plot for two particles' trajectories with interactions in the (z,vz) plane
    axs[0, 1].plot(z1_w,vz1_w, color='salmon', label='Particle 1')
    axs[0, 1].plot(z2_w,vz2_w, color='red', label='Particle 2')
    axs[0, 1].set_title('... in the (z,$v\_{z}$) plane, with interactions')
    axs[0, 1].xlabel('z [\u03bcm]')
    axs[0, 1].ylabel('v\_{z}[\u03bcm/\u03bcs]')
    axs[0, 1].legend()

    # Plot for two particles' trajectories without interactions in the (x,vx) plane
    axs[1, 0].plot(x1_wo, vx1_wo, color='teal', label='Particle 1')
    axs[1, 0].plot(x2_wo, vx2_wo, color='mediumturquoise', label='Particle 2')
    axs[1, 0].set_title('... in the (x,$v\_{x}$) plane, without interactions')
    axs[1, 0].xlabel('x [\u03bcm]')
    axs[1, 0].ylabel('$v\_{x}$ [\u03bcm/\u03bcms]')
    axs[1, 0].legend()
     
    # Plot for two particles' trajectories without interactions in the (z,vz) plane
    axs[1, 1].plot(z1_wo, vz1_wo, color='darkmagneta',label='Particle 1')
    axs[1, 1].plot(z2_wo, vz2_wo, color='orchid', label='Particle 2')
    axs[1, 1].set_title('... in the (x,$v\_{x}$) plane, without interactions')
    axs[1, 1].xlabel('z [\u03bcm]')
    axs[1, 1].ylabel('v\_{z}[\u03bcm/\u03bcs]')
    axs[1, 1].legend()

    # Saving the plots as an svg file and opening it
    plt.savefig("two_particle_traj.svg", format="svg") 
    plt.show()


    return 0