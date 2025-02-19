import numpy as np
import matplotlib.pyplot as plt
from scipy import constants

from custom_functions.read import read_observables

t_au, Pt, Et_au, _ = read_observables("./", "new")
t_fs = t_au*constants.physical_constants["atomic unit of time"][0]*1.e15

fig, ax = plt.subplots(2)
for ip,p in enumerate(Pt):
#    ax[0].plot(t_fs, np.log(p), label=str(ip))
    ax[0].plot(t_fs, p, label=str(ip))
ax[0].legend()
#ax[0].set_xticks([])
#ax[0].set_ylim(0.00087, 0.00089)
#ax[0].set_ylim(-1.e-05, 1.e-05)
ax[1].plot(t_fs, Et_au[0])
ax[1].plot(t_fs, Et_au[1])
ax[1].plot(t_fs, Et_au[2])
ax[1].set_xlabel("Time (fs)")
plt.show()
plt.savefig("Population.png", dpi=200)
