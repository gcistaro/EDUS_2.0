import numpy as np
import matplotlib.pyplot as plt
from scipy import constants

from custom_functions.read import read_observables

import argparse, sys

#define arguments
parser=argparse.ArgumentParser()

parser.add_argument("--folder", default="./Output/", help="Folder where the .txt files are contained")
parser.add_argument("--version", default="new", help="Version of the EDUS code that you used")
args=parser.parse_args()

#print in output some infos
print("folder       :    ", args.folder)
print("version      :    ", args.version)
t_au, Pt, Et_au, Vt_au = read_observables(args.folder, args.version)
t_fs = t_au*constants.physical_constants["atomic unit of time"][0]*1.e15

#fig, ax = plt.subplots(3, figsize=(10,10))
fig, ax = plt.subplots(3)
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
ax[2].plot(t_fs, np.real(Vt_au[0]))
ax[2].plot(t_fs, np.real(Vt_au[1]))
ax[2].plot(t_fs, np.real(Vt_au[2]))
ax[2].set_xlabel("Time (fs)")
#plt.savefig("Population.pdf")
plt.savefig("Population.png", dpi=200)
plt.show()
