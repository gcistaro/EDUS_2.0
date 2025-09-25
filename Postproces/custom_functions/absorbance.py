
import numpy as np
from custom_functions.fouriertransform import FourierTransform
from scipy import constants
import matplotlib.pyplot as plt

def absorbance(t_au, Vt_au, Et_au, limits, smearing = 0., plot = False, data = False):
    # print("Fourier transforming")
    freq_eV, Jw_au = FourierTransform(t_au, Vt_au, smearing)
    _, Ew_au       = FourierTransform(t_au, Et_au)
    # print("done")


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
        # index = np.logical_and(np.logical_and( freq_eV > limits[0], freq_eV < limits[1] ), np.abs(Ew_au) > 1e-14)
        # print(index)
        # index = np.where(index)[1]
        # print(index)

    Jw_au = Jw_au[:,index]
    Ew_au = Ew_au[:,index]
    freq_eV = freq_eV[index]
    
    freq_au = freq_eV*constants.physical_constants["hartree-electron volt relationship"][0]

    alpha = constants.physical_constants["fine-structure constant"][0]
    
    dw_au = (1j * Jw_au) /freq_au #Compute the dipole spectrum form d(w) = -z(w) and v = dz/dt => v(w) =  iwz(w) and d(w) = iv(w)/w

    Sw_au  = -(2*np.imag( np.sum( [ np.conj(Ew_au[ix])*dw_au[ix] for ix in [0,1,2] ], axis=0) ) )        #response function

    Iw_au = 2*np.sum( [np.abs(Ew_au[ix])*np.abs(Ew_au[ix]) for ix in [0,1,2] ], axis=0)/(4.*np.pi*alpha*freq_au)

    ######## ONLY USE TO PREVENT DIVISION BY 0 ##########
    for freq in range(len(Iw_au)):                    #
                if np.abs(Iw_au[freq]) < 1e-14:       #
                    Iw_au[freq] = 0 # set to true 0   #
    #####################################################

    Absorbance =np.divide(Sw_au,Iw_au, np.zeros_like(Sw_au, dtype = float), where=Iw_au!=0) # If divide by 0, Absorbance is set to 0.


    if plot:
        plt.plot(freq_eV, np.abs(dw_au[0]), label = r"Along x")
        plt.plot(freq_eV, np.abs(dw_au[1]), label = r"Along y")
        plt.plot(freq_eV, np.abs(dw_au[2]), label = r"Along z")
        plt.legend()
        plt.show()
    
    if data:
        print("Ew[1] min and max: ", np.min(np.abs(Ew_au[1])), np.max(np.abs(Ew_au[1])))
        print("Ew[2] min and max: ", np.min(np.abs(Ew_au[2])), np.max(np.abs(Ew_au[2])))
        print("Ew min and max: ", np.min(np.abs(Ew_au)), np.max(np.abs(Ew_au)))
        print("Vw min and max: ", np.min(np.abs(Jw_au)), np.max(np.abs(Jw_au)))
        print("rw min and max: ", np.min(np.abs(dw_au)), np.max(np.abs(dw_au)))
        print("Sw min and max: ", np.min(np.abs(Sw_au)), np.max(np.abs(Sw_au)))
        print("Iw min and max: ", np.min(np.abs(Iw_au)), np.max(np.abs(Iw_au)))
        print("Absorbance min and max: ", np.min(np.abs(Absorbance)), np.max(np.abs(Absorbance)))
        # used for debugging :
        # print(Iw_au)
        # print(Absorbance)
    

    


    return freq_eV, Absorbance



