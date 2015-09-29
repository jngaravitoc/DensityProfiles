import matplotlib.pyplot as plt
import numpy as np


data = np.loadtxt('orbits.txt')

t = data[:,0]
x = data[:,1]
y = data[:,2]
z = data[:,3]
vx = data[:,4]
vy = data[:,5]
vz = data[:,6]


r = np.sqrt(x**2 + y**2 + z**2)

plt.plot(t, r, lw = 2, alpha=0.7)
plt.show()

