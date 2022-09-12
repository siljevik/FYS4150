import matplotlib.pyplot as plt

#Creating empty lists
x = []
u = []
# Opening the .txt-file created by the c++ code
with open('x_u.txt') as f:
    # Reading the file
    data = f.read()
    # Defining each line
    lines = data.split('\n')
    n = len(lines) # length of lines list
    i = 0
    m = n-2 # The last two lines in the lines list are useless
    # So we only look at the lines that are not the two last ones
    for k in range(m):
        line = lines[k] # k as in line number
        list_x_u = line.split()
        x_v = list_x_u[0]
        u_v = list_x_u[1]
        x.append(float(x_v))
        u.append(float(u_v))      
    list_last_x_u = lines[n-2].split()
    x_last_v = list_last_x_u[0]
    u_last_v = list_last_x_u[1]
    x.append(float(x_last_v))
    u.append(0) # Somehow we got the wrong number here, so we set it to the correct, 0
    
    # Then we can plot it
    plt.plot(x,u)
    # Labeling plot, x and y axes
    plt.title('Plot of Poisson equation')
    plt.xlabel('x')
    plt.ylabel('u(x)')
    plt.show()