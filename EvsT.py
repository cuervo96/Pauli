import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("EvsT_rho=0.08.txt")
data1 = np.loadtxt("EvsT_rho=0.10.txt")
data2 = np.loadtxt("EvsT_rho=0.12.txt")
data3 = np.loadtxt("EvsT_rho=0.15.txt")
data4 = np.loadtxt("EvsT_rho=0.16.txt")
data5 = np.loadtxt("EvsT_rho=0.17.txt")
data6 = np.loadtxt("EvsT_rho=0.18.txt")
j = data[:,0]
E = data[:,1]
j1 = data1[:,0]
E1 = data1[:,1]
E2 = data2[:,1]
E3 = data3[:,1]
E4 = data4[:,1]
E5 = data5[:,1]
E6 = data6[:,1]
T = np.zeros(len(j))
T1 = np.zeros(len(j))
T1[0] = 4
for i in range(len(j)):
    T[i] = 4 - (4-0.05)* j[i] / 1000
for i in range(len(j1) - 1):
    T1[i + 1] = 4- (4-0.05)* j1[i] / 100
plt.figure()
plt.plot(T1[1:len(j1)],E[1:len(E)]/512)
plt.plot(T1[1:len(j1)],E1[1:len(E)]/512)
plt.plot(T1[1:len(j1)],E2[1:len(E)]/512)
plt.plot(T1[1:len(j1)],E3[1:len(E)]/512)
plt.plot(T1[1:len(j1)],E4[1:len(E)]/512)
plt.plot(T1[1:len(j1)],E5[1:len(E)]/512)
plt.plot(T1[1:len(j1)],E6[1:len(E)]/512)
plt.grid()
plt.title("Curva Calórica")
plt.xlabel("T (MeV)")
plt.ylabel("E (MeV)")
plt.legend(["ρ = 0.08","ρ = 0.10","ρ = 0.12","ρ = 0.15","ρ = 0.16","ρ = 0.17", "ρ = 0.18"])
plt.savefig("Pauli.png")
plt.show()
