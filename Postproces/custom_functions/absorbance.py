
import numpy as np
from custom_functions.fouriertransform import FourierTransform
from scipy import constants
import matplotlib.pyplot as plt

def absorbance(t_au, Vt_au, Et_au, limits):
    print("Fourier transforming")
    freq_eV, Jw_au = FourierTransform(t_au, Vt_au, True)
    _, Ew_au       = FourierTransform(t_au, Et_au, False)
    print("done")


    #cut a window
    index = np.where(np.logical_and( freq_eV > limits[0], freq_eV < limits[1] ))[0]
    Jw_au = Jw_au[:,index]
    Ew_au = Ew_au[:,index]
    freq_eV = freq_eV[index]
    freq_au = freq_eV*constants.physical_constants["hartree-electron volt relationship"][0]

    freq_au = freq_eV*constants.physical_constants["hartree-electron volt relationship"][0]
    alpha = constants.physical_constants["fine-structure constant"][0]
    
    rw_au = 1j*freq_au*Jw_au
    Sw_au  = -(2*np.imag( np.sum( [ np.conj(Ew_au[ix])*rw_au[ix] for ix in [0,1,2] ], axis=0) ) )        #response function
    Iw_au = 2*np.sum( [np.abs(Ew_au[ix])*np.abs(Ew_au[ix]) for ix in [0,1,2] ], axis=0)/(4.*np.pi*alpha*freq_au)
    plt.plot(freq_au, rw_au[0])
    plt.show()
    
    print("Ew[1] min and max: ", np.min(np.abs(Ew_au[1])), np.max(np.abs(Ew_au[1])))
    print("Ew[2] min and max: ", np.min(np.abs(Ew_au[2])), np.max(np.abs(Ew_au[2])))
    print("Ew min and max: ", np.min(np.abs(Ew_au)), np.max(np.abs(Ew_au)))
    print("Vw min and max: ", np.min(np.abs(Jw_au)), np.max(np.abs(Jw_au)))
    print("rw min and max: ", np.min(np.abs(rw_au)), np.max(np.abs(rw_au)))
    print("Sw min and max: ", np.min(np.abs(Sw_au)), np.max(np.abs(Sw_au)))
    print("Iw min and max: ", np.min(np.abs(Iw_au)), np.max(np.abs(Iw_au)))
    Absorbance =Sw_au/Iw_au

    return freq_eV, Absorbance
