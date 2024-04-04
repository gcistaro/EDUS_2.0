import numpy as np
import matplotlib.pyplot as plt



Population = np.loadtxt("Population.txt", unpack=True)
t = range(Population.shape[-1])
Laser = np.loadtxt("Laser.txt", unpack=True)
Velocity = np.loadtxt("Position.txt", unpack=True)


fig, ax = plt.subplots(2)
for p in Population:
    ax[0].plot(t, p)
ax[0].plot(t,np.zeros(Population.shape[-1]))
#ax[0].set_ylim(-1.e-05, 1.e-05)
ax[1].plot(t, Laser[0])
ax[1].plot(t, Laser[1])
ax[1].plot(t, Laser[2])

plt.show()
plt.savefig("Population.png", dpi=200)
