import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("test_metropolis.txt")

plt.figure()
plt.plot(data[0:len(data)-1])
plt.grid()
plt.xlabel("Pasos/100")
plt.ylabel("Energia (MeV)")
plt.legend(["T = 2 Mev"])
plt.show()
