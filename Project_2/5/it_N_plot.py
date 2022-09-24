## plot number of terations against size of the matrix
## some manual labor required to get the values because we didn't make
## a loop in our C++ code:)

import matplotlib.pyplot as plt
  
x = []
y = []
for line in open('iterations_dense.txt', 'r'):
    lines = [i for i in line.split()]
    x.append(lines[0])
    y.append(int(lines[1]))
     
plt.title("Required iterations- for a dense matrix")
plt.xlabel('Matrix size [N]')
plt.ylabel('Iterations')
plt.yticks(y)
plt.plot(x, y, marker = 'o', c = 'g')
plt.savefig('line_plot_dense.pdf')   
plt.show()
