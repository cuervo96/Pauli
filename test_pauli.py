import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("test_metropolis.txt")

plt.figure()
plt.plot(data[0:len(data)-1], '.')
plt.grid()
plt.show()
