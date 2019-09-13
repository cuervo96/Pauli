import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("test_Tablas.txt")
#pauli = data[:,0]
LJ = data

plt.figure()
plt.plot(LJ, '.')
plt.grid()
plt.xlim([0,len(LJ)])
plt.ylim([-10,10])
plt.show()

