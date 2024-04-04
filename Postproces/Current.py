#from scipy.fft import fft, fftfreq, fftshift
from scipy import constants
import numpy as np
from numpy import fft
import matplotlib.pyplot as plt

#read data
J = np.loadtxt("Position.txt", unpack=True)
Jx_t = J[0]
Jy_t = J[2]
Jz_t = J[4]

length_fft=20000
Jy_w = fft.fftshift(fft.fft(Jy_t, n=length_fft))
freq = fft.fftshift(fft.fftfreq(length_fft, d=0.4*constants.physical_constants["atomic unit of time"][0]))
freq = freq/constants.physical_constants["electron volt-hertz relationship"][0]
print(freq)
plt.scatter(freq, np.abs(Jy_w))#.real, freq, Jy_w.imag)
plt.xlim(-10,10)
#plt.ylim(-0.01,0.01)
plt.show()
plt.savefig("velocity.png", dpi=200)
