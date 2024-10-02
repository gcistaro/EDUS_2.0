#from scipy.fft import fft, fftfreq, fftshift
from scipy import constants
from scipy.constants import e, hbar
import numpy as np
from numpy import fft
import matplotlib.pyplot as plt
import sys

#FourierTransform function 
def FourierTransform(t, Jt, cut):
    Delta_t_s = ((t[1]-t[0])*constants.physical_constants["atomic unit of time"][0])
    T_s = (t[-1]-t[0])*constants.physical_constants["atomic unit of time"][0]
    Delta_freq_eV_wanted = 0.01
    T_s_wanted = 1./Delta_freq_eV_wanted/constants.physical_constants["electron volt-hertz relationship"][0]
    length_fft = int( T_s_wanted/Delta_t_s )
    print(length_fft)
    newt = np.linspace(t[0],T_s_wanted,length_fft)
    
    #alpha to use to get a decent one body part:
    #alpha = -(t[-1]-t[0])/np.log(1.e-02)

    alpha = -(t[-1]-t[0])/np.log(1.e-01)
    print(np.exp(-t[0]/alpha), np.exp(-t[-1]/alpha))
    if(cut) :
        window = np.exp(-t/alpha)#np.blackman(length_fft) 
        #window=1
        #window = np.blackman(t.shape[0]) 
    else:
        window = 1. 
    print(window)
    plt.plot(t, Jt[0], label="nowindow")
    plt.plot(t, window*Jt[0], label="window")
    plt.legend()
    plt.show()
    Jw = [fft.fftshift(fft.ifft(window*Jt[ix], n=length_fft)) for ix in np.array([0,1,2])]
    freq = fft.fftshift(fft.fftfreq(length_fft, d=(t[1]-t[0])*constants.physical_constants["atomic unit of time"][0]))
    freq = freq/constants.physical_constants["electron volt-hertz relationship"][0]
    return freq, Jw


#read data
t=np.loadtxt("Time.txt")
Et = np.loadtxt("Laser.txt", unpack=True)
Vt = np.loadtxt("Velocity.txt", unpack=True)
Vt = np.array([Vt[0], Vt[2], Vt[4]])
freq, Ew = FourierTransform(t, Et, False)
freq, Vw = FourierTransform(t, Vt, True)
rw = 1j*freq*Vw

Sw  = -(2*np.imag(np.conj(Ew[0])*rw[0]))
Iw = 2*np.abs(Ew[0])*np.abs(Ew[0])/freq
print(np.max(Sw))
Absorption = Sw/Iw 
index = np.where(np.logical_and( freq > 4.4, freq<12))#7.25) ) 
Ew[0] = Ew[0][index]
Ew[1] = Ew[1][index]
Ew[2] = Ew[2][index]
Sw = Sw[index]
freq = freq[index]
Absorption = Absorption[index]
plt.plot(freq, np.abs(Ew[0][:])*np.max(Sw)/np.max(np.abs(Ew[0][:])), label="E"+str(0))
plt.plot(freq, Sw, label="Sw" )
plt.plot(freq, Absorption*np.max(Sw)/np.max(np.abs(Absorption)), label="Abs(w)")
plt.axvline(x=7.25)
plt.ylim(-0.1*np.max(Sw),np.max(Sw))
plt.legend()
plt.show()
#plt.savefig("velocity.png", dpi=200)
