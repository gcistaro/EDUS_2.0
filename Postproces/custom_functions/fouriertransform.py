import numpy as np
from scipy import constants
from numpy import fft



#FourierTransform function 
def FourierTransform(t_au, Jt_au, cut):
    Delta_t_s = ((t_au[1]-t_au[0])*constants.physical_constants["atomic unit of time"][0])
    Total_T_s = (t_au[-1]-t_au[0])*constants.physical_constants["atomic unit of time"][0]
    Delta_freq_eV_wanted = 0.001
    Total_T_s_wanted = 1./Delta_freq_eV_wanted/constants.physical_constants["electron volt-hertz relationship"][0]
    length_fft = int( Total_T_s_wanted/Delta_t_s )

    print(Total_T_s_wanted,Delta_t_s, Total_T_s,Delta_freq_eV_wanted)
    t_au_wanted = np.linspace(t_au[0],Total_T_s_wanted,length_fft)
    
    #alpha to use to get a decay:
    alpha = -(t_au[-1]-t_au[0])/np.log(1.e-04)

    if(cut) :
        window = np.exp(-t_au/alpha)
    else:
        window = 1. 

    Jw_au = -np.array([fft.fftshift(fft.ifft(window*Jt_au[ix], n=length_fft)) for ix in np.array([0,1,2])])
    freq_Hz = fft.fftshift(fft.fftfreq(length_fft, d=(t_au[1]-t_au[0])*constants.physical_constants["atomic unit of time"][0]))
    freq_eV = freq_Hz/constants.physical_constants["electron volt-hertz relationship"][0]
    return freq_eV, Jw_au
