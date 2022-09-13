import numpy as np
import math as m
import matplotlib.pyplot as plt

n = 10
h = 1/n

###### b
bs = []
b_n = 2

for i in range(n+1): # i går opp til, men ikke til og med n+1
    bs.append(b_n)
    # lager ny b
    b = 2 - (1/b_n) 
    # oppdaterer for ny b
    b_n = b

##### g
gs = []
x = 0

for i in range(n+1):
    x += h*i 
    g = (h**2)*100*m.e**(-10*x)
    gs.append(g)


##### v
v=[]
v_n = 0 # Because of boundary conditions and that we do it backwards
v.append(v_n)
for i in range(n,-1,-1):    
    vn = (gs[i]+v_n)/(bs[i])
    v.append(vn)
    v_n = vn # Updating the v
v.reverse() # since we did it backwards


xs = []
x=1 # Starting backwards here too
xs.append(x)
for i in range(n,-1,-1):
    x -= h
    xs.append(x)
xs.reverse() # since we did it backwards

# Then we can plot it
plt.plot(xs,v)
# Labeling plot, x and y axes
plt.xlabel('x')
plt.ylabel('v(x)')
# Saving the figure as a pdf-file
plt.savefig("plot_7_10.pdf")
plt.show()