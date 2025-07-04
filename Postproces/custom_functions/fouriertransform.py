import numpy as np
from scipy import constants
from numpy import fft



def FourierTransform(t_au, Jt_au, smearing_eV):
    """
    This function calculates the fourier transform of Jt_au given the time grid t_au. 
    The sign convention is: 
    J(w) = \int dt J(w)e^{+iwt} 
    It is also possible to define a smearing for the function in time. 
    This means that the function that will be fourier transformed is:
    J(t) = J(t)*exp(-t*smearing)
    The default units for the smearing are electronvolts 
    """
    Delta_t_s = ((t_au[1]-t_au[0])*constants.physical_constants["atomic unit of time"][0])
    Total_T_s = (t_au[-1]-t_au[0])*constants.physical_constants["atomic unit of time"][0]
    Delta_freq_eV_wanted = 0.001
    Total_T_s_wanted = 1./Delta_freq_eV_wanted/constants.physical_constants["electron volt-hertz relationship"][0]
    length_fft = int( Total_T_s_wanted/Delta_t_s )

    #print(Total_T_s_wanted,Delta_t_s, Total_T_s,Delta_freq_eV_wanted)
    t_au_wanted = np.linspace(t_au[0],Total_T_s_wanted,length_fft)
    
    smearing_au = smearing_eV*constants.physical_constants["electron volt-hartree relationship"][0]

    window = np.exp(-t_au*smearing_au)

    Jw_au = -np.array([fft.fftshift(fft.ifft(window*Jt_au[ix], n=length_fft)) for ix in np.array([0,1,2])])
    freq_Hz = fft.fftshift(fft.fftfreq(length_fft, d=(t_au[1]-t_au[0])*constants.physical_constants["atomic unit of time"][0]))
    freq_eV = freq_Hz/constants.physical_constants["electron volt-hertz relationship"][0]
    return freq_eV, Jw_au
