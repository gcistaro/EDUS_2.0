from custom_functions.read import *
from custom_functions.fouriertransform import *
import matplotlib.pyplot as plt


def get_emission(datadir):
    '''Calculates emission spectra. 
    
    Parameters
    ----------
    datadir: str
        Path to data folder.
    '''
     
    t_au, _, _, Vt_au = read_observables(datadir, "new") # We only need velocity

    freq_eV, Vw_au = FourierTransform(t_au, Vt_au) # Do the FT 
    freq_au = freq_eV*constants.physical_constants["hartree-electron volt relationship"][0]

    Aw_au = Vw_au * 1j*freq_au 

    return freq_eV, Aw_au

def plotHHG(datadir, fund,start = 0, stop=10, axes=[]):
    """Plots the emission spectra, in semi-log scale, from the Fourier transform of the velocity. 

    Parameters
    ----------
    datadir: str
        Path to data folder.
    fund: float
        "Fundamental" frqeuncy of the laser, in eV. Used to set the x axis.
    start: int, optional
        First harmonic order plotted.
    stop: int, optional
        Highest harmonic order plotted.
    axes: list of str, optional
        Direction of the acceleration to plot. If left empty, will plot the whole norm of the acceleration.
    
    Returns
    -------
    fig: matplotlib.pyplot.Figure
        The spectra in semilog scale.
    
    Raises
    ------
    ValueError
        If the axes are not provided properly, i.e. different than x, y and z or empty.

    Notes
    -----
    The emission spectra is calculated from :math:`S(\\omega) \\propto \\omega^2 \\Vert v(\\omega) \\Vert ^2` where :math:`v(\\omega)` is the velocity[1]_.

    .. [1] Lun Yue and Mette B. Gaarde, *Introduction to theory of high-harmonic generation in solids: tutorial* J. Opt. Soc. Am. B **39** (2022)
    
    """
    ###### Check inputs ######
    if ('x' not in axes) and ('y' not in axes) and ('z' not in axes) and axes!=[]:
            raise ValueError(f"Axis are poorly defined")
    

    freq_eV, Aw_au = get_emission(datadir) #Compute the emission

    Aw_au_x = Aw_au[0,:]
    Aw_au_y = Aw_au[1,:]
    Aw_au_z = Aw_au[2,:]

    Aw = np.linalg.norm(Aw_au, axis=0) #norm over all axis

    ###### plotting ######
    fig, ax = plt.subplots()
    if axes==[]:
        ax.semilogy(freq_eV/fund, Aw**2)
    
    for axe in axes:
        ax.semilogy(freq_eV/fund, np.abs(vars()["Aw_au_" + axe])**2, label = r"Along " + axe)

    ax.set_xticks(np.linspace(start,stop, stop +1-start, dtype = int)) # get all the orders
    ax.set_xlim(start,stop)
    ax.set_xlabel(r"Harmonic order")
    ax.set_ylabel(r"Intensity (arb. units, log scale)")
    ax.legend()

    return fig