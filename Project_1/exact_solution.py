import matplotlib.pyplot as plt
import numpy as np
import math

def exact_solution():
    n = 100
    x = np.arange(0, n, 1) # (start, stop, stepsize)
    y = 1-(1-math.e^(-10))*x-math.e^(-10*x)
    plt.plot(y)

    plt.show()