
import numpy as np
from custom_functions.fouriertransform import FourierTransform
from scipy import constants

def get_absorbance(t_au, Vt_au, Et_au, limits=[], smearing=0.):
    """
    This function calculates the absorbance in frequency domain given as input the time grid (WARNING, need to be equally spaced),
    The velocity over the time grid, the electric field using for the time propagation. 
    You can also define limits where you want the absorbance in units of eV (is frequency) and the smearing. 
    For more details about the smearing, please check the function FourierTransform.
    """

    freq_eV, Vw_au = FourierTransform(t_au, Vt_au, smearing)
    _, Ew_au       = FourierTransform(t_au, Et_au, 0.)

    #cut a window
    index = []
    if limits == []:
        #take only frequencies for which the EF is at least 5% of its maximum
        #and frequency is positive
        Ew_norm = np.linalg.norm(Ew_au[:,:], axis=0)
        index = np.where( np.logical_and(
                                            Ew_norm >= 0.05*np.max(Ew_norm),
                                            freq_eV > 0.0 )
                        )[0]
    else:
        #take only values in a window of frequency that you can define in the input
        index = np.where(np.logical_and( freq_eV > limits[0], freq_eV < limits[1] ))[0]

    #define everything only in the window
    Vw_au   = Vw_au[:,index]
    Ew_au   = Ew_au[:,index]
    freq_eV = freq_eV[index]
    freq_au = freq_eV*constants.physical_constants["hartree-electron volt relationship"][0]

    freq_au = freq_eV*constants.physical_constants["hartree-electron volt relationship"][0]
    alpha = constants.physical_constants["fine-structure constant"][0]

    dw_au = Vw_au/(1j*freq_au)

    #response function
    Sw_au  = (2*np.imag( np.sum( [ np.conj(Ew_au[ix])*dw_au[ix] for ix in [0,1,2] ], axis=0) ) )        #response function
    
    #spectral function of the electric field
    Iw_au = 2*np.sum( [np.abs(Ew_au[ix])*np.abs(Ew_au[ix]) for ix in [0,1,2] ], axis=0)/(4.*np.pi*alpha*freq_au) # not sure about these constants
    Absorbance =Sw_au/Iw_au

    return freq_eV, Absorbance
