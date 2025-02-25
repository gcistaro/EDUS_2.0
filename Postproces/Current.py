#from scipy.fft import fft, fftfreq, fftshift
from scipy import constants
import numpy as np
from numpy import fft
import matplotlib.pyplot as plt
import sys

from custom_functions.read import read_observables
from custom_functions.absorbance import absorbance



#read data
folder = sys.argv[1] 

if len(sys.argv) > 2:
    version = sys.argv[2]
else:
    version = "new"

t_au, _, Et_au, Vt_au = read_observables(folder, version)

plt.plot(t_au, Vt_au[0])
plt.show()
freq_eV, Absorbance = absorbance(t_au, Vt_au, Et_au,[1,7])


plt.plot(freq_eV, Absorbance, label="Abs(w)")
#plt.axvline(x=7.25)
plt.legend()
plt.show()
np.savetxt("absorbance.txt", np.array(np.transpose([freq_eV, Absorbance])))
plt.savefig("Absorbance.png", dpi=200)
