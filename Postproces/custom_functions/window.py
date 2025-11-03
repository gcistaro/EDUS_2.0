import numpy as np

def CutWindow(freq_eV, Ew_norm, limits=None):
    #cut a window
    print("Defining window in energy (eV)...", end=" ", flush=True)
    index = []
    if limits == None:
        #take only frequencies for which the EF is at least 5% of its maximum
        #and frequency is positive
        index = np.where( np.logical_and(
                                            Ew_norm >= 0.05*np.max(Ew_norm),
                                            freq_eV > 0.0 )
                        )[0]
    else:
        #take only values in a window of frequency that you can define in the input
        index = np.where(np.logical_and( freq_eV > limits[0], freq_eV < limits[1] ))[0]
    print("done.")
    return index

