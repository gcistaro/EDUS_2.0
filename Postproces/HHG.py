#from scipy.fft import fft, fftfreq, fftshift
from scipy import constants
import numpy as np
from numpy import fft
import matplotlib.pyplot as plt
import sys

from custom_functions.read import read_observables
from custom_functions.fouriertransform import FourierTransform



#read data
folder = sys.argv[1] 

if len(sys.argv) > 2:
    version = sys.argv[2]
else:
    version = "new"

t_au, _, _, Vt_au = read_observables(folder, version)

plt.plot(t_au, Vt_au[0])
plt.show()
freq_eV, Vw = FourierTransform(t_au, Vt_au, True)

#rw = Vt_au/(1j*freq_eV)
Vw = np.linalg.norm(Vw, axis=0)
print("Vw shape:" , Vw.shape)


plt.plot(freq_eV, np.log10(Vw), label="$log_{10}(r(\omega))$")
#plt.axvline(x=7.25)
plt.legend()
plt.show()
np.savetxt("absorbance.txt", np.array(np.transpose([freq_eV, Absorbance])))
plt.savefig("Absorbance.png", dpi=200)
