from custom_functions.read import *
from custom_functions.fouriertransform import *
import matplotlib.pyplot as plt
from custom_functions.absorbance import absorbance
from scipy.signal import find_peaks
from custom_functions.Exp_manager import * 

def assertAxes(axes):
    if ('x' not in axes) and ('y' not in axes) and ('z' not in axes) and axes!=[]:
            raise ValueError(f"Axis are poorly defined.")

def plotHHG(exp, start=0, stop=10, axes = [], smearing=0., **kwargs):
    assertAxes(axes)
    
    exp.gen_HHG_data(smearing=smearing)

    x_harmonic_axis = exp.frequ_axis_eV/exp.HHG_fund_freq_eV

    fig, ax = plt.subplots()

    if axes == []:
        ax.plot(x_harmonic_axis, exp.HHG_acceleration_freq_au**2, **kwargs)
    
    if 'label' in kwargs.keys():
        del kwargs['label'] #LAbel are imposed in the axis plotting.

    for axe in axes:
        ax.plot(x_harmonic_axis, np.abs(vars(exp)["HHG_acceleration_freq_au_"+ axe] )**2, label = r"Along " + axe, **kwargs)

    ax.set_yscale('log')
    ax.set_xticks(np.linspace(start,stop, stop +1-start, dtype = int)) # get all the orders
    ax.set_xlim(start,stop)
    ax.set_xlabel(r"Harmonic order")
    ax.set_ylabel(r"Intensity (arb. units, log scale)")
    ax.legend()

    return fig

def plotHHGPhase(exp, start=0, stop=10, smearing=0., axis = 0,  **kwargs):
    # assertAxes(axes)
    
    exp.gen_HHG_peaks(smearing=smearing)

    x_harmonic_axis = exp.frequ_axis_eV/exp.HHG_fund_freq_eV
    Vw_au = exp.velocity_FT_dict[f'{smearing}']

    V = np.dot(axis, Vw_au)
    delta_phi = []

    phase = np.angle(V)
    indices = exp.HHG_harm_idx
    for i in range(len(indices)-1):
        delta_phi.append(phase[indices[i+1]]-phase[indices[i]])

    fig, ax = plt.subplots()

    # if axes == []:
    ax.plot(x_harmonic_axis, phase, **kwargs)
    
    if 'label' in kwargs.keys():
        del kwargs['label'] #LAbel are imposed in the axis plotting.

    # for axe in axes:
    #     ax.plot(x_harmonic_axis, np.abs(vars(exp)["HHG_acceleration_freq_au_"+ axe] )**2, label = r"Along " + axe, **kwargs)

    # ax.set_yscale('log')
    ax.set_xticks(np.linspace(start,stop, stop +1-start, dtype = int)) # get all the orders
    ax.set_xlim(start,stop)
    ax.set_xlabel(r"Harmonic order")
    ax.set_ylabel(r"Phase (radians)")
    ax.legend()

    return fig, delta_phi


def plotHHGPeaks(exp, start=0, max_harmonic = 20, smearing= 0., cutoff= 1.5e-1, normed = False, **kwargs):

    exp.gen_HHG_peaks(max_harmonic=max_harmonic, smearing=smearing, cutoff = cutoff)

    x_harmonic_axis = exp.frequ_axis_eV/exp.HHG_fund_freq_eV
    harm_idx = exp.HHG_harm_idx

    fig, ax = plt.subplots()

    if normed:
        y_axis = exp.HHG_normed_pics_height
    else:
        y_axis = exp.HHG_pics_height
    ax.scatter(x_harmonic_axis[harm_idx], y_axis, **kwargs)

    ax.set_yscale('log')

    ax.set_xticks(np.linspace(start,max_harmonic, max_harmonic-start+1, dtype = int)) # get all the orders
    ax.set_xlim(start,max_harmonic)
    ax.set_xlabel(r"Harmonic order")
    ax.set_ylabel(r"Intensity (arb. units, log scale)")
    ax.legend()

    return fig

def get_peak_intesity(exp, order, smearing= 0., cutoff= 1.5e-1, normed = False):

    exp.gen_HHG_peaks(max_harmonic = order+1, smearing = smearing, cutoff = cutoff)

    if len(exp.HHG_pics_height) < order and order > 5:
        print("Peak not found.")
        return -1
    
    x = exp.frequ_axis_eV/exp.HHG_fund_freq_eV

    if normed:
        peaks = exp.HHG_normed_pics_height
    else: 
        peaks = exp.HHG_pics_height

    idx = np.argmin(np.abs(x[exp.HHG_harm_idx] -  order))

    return peaks[idx]

def plotPopulation(exp, **kwargs):
     
    t_au = exp.time_au
    Pt = exp.population_time
    nb_bands = exp.filled_bands

    occupied_bands = [] #initially occupied, shows the number of holes
    empty_bands = [] #initially empty, shows the number of electrons

    for i, band in enumerate(Pt):
        if np.sum(band) >= 1e-6:
            if i <= nb_bands-1:
                occupied_bands.append(1-band)
            else: 
                empty_bands.append(band)
    

    fig, ax = plt.subplots(len(occupied_bands)+len(empty_bands), sharex=True)
    i = 0
    for band in empty_bands+occupied_bands:
        ax[i].plot(t_au, band, label = f'Band n° %s' % i, color = f'C{i}', **kwargs)
        # ax[i].text(0.10, 0.80, f'Band n° %s' % i, horizontalalignment='center', transform=ax[i].transAxes, fontsize=12)
        i+= 1

    return fig

def plotVelocity(exp, space='time', axes=[], time_window = (0,0), freq_lims= (0,15), smearing = 0.):
    """ Plots the velecity in fuction of space or frequency alongsides the given directions.

    Parameters
    ----------
    exp : instance of Experience
        EDUS experience.
    space : str, optional
        Space of plotting, can be 'time', 'frequency' or 'both'.
    axes : list of str, optional
        Direction of current to plot. If left empty, will plot the whole norm of the current.
    time_window : tuple of floats, optional
        Time interval of the simulation to consider in au. If none provided, Fourier transform is computed on the whole simulation. 

    Returns
    -------
    fig: matplotlib.pyplot.Figure
        The absorption spectra.

    Raises
    -------
    ValueError
        If the time_window provided is not in proper order.
    ValueError
        If the axes are not provided properly, i.e. different than x, y and z or empty.
    ValueError
        If the space provided is not 'time', 'frequency' or 'both'.
    """

    ###### Check inputs ######
    if time_window[0] > time_window[1]:
        raise ValueError(f"You cannot reverse time !")
    
    if ('x' not in axes) and ('y' not in axes) and ('z' not in axes) and axes!=[]:
            raise ValueError(f"Axis are poorly defined")
    
    if space not in {'time', 'frequency', 'both'}:
         raise ValueError(f'Plotting space must be one of %r.' % {'time', 'frequency', 'both'})
    
    exp.gen_observables(smearing = smearing)

    t_au = exp.time_au #read output
    V_au = exp.velocity_time_au
    

    ###### Select the proper plotting space ######
    time, freq = False, False

    if space == 'time':
        time = True
        fig,ax = plt.subplots()
        ax=[ax]
    elif space == 'frequency':
        freq = True
        fig,ax = plt.subplots()
        ax=[ax]
    elif space == 'both':
        time, freq = True, True
        fig,ax = plt.subplots(2)

    ###### Sort the data accordingly ######
    if time:
        V_au_x = V_au[0,:]
        V_au_y = V_au[1,:]
        V_au_z = V_au[2,:]

    if freq:
        freq_axis = exp.frequ_axis_eV
        Vw_au = exp.velocity_FT_dict[f'{smearing}']

        Vw_au_x = Vw_au[0,:]
        Vw_au_y = Vw_au[1,:]
        Vw_au_z = Vw_au[2,:]
         
    ###### Plotting ######

    if axes == []:
        i = 0
        if time:
            ax[i].plot(t_au, np.linalg.norm(V_au, axis = 0))
            ax[i].set_xlabel(r"Time (in a.u.)")
            ax[i].set_ylabel(r"Velocity (in a.u.)")
            i+= 1
        if freq:
            ax[i].plot(freq_axis, np.linalg.norm(Vw_au, axis= 0))
            ax[i].set_xlim(freq_lims)
            # ax[i].
        
    for axe in axes:
        i=0
        if time:
            ax[i].plot(t_au, vars()['V_au_' + axe], label = r"Along "+ axe)
            i+=1
        if freq:
            ax[i].plot(freq_axis, vars()['Vw_au_'+axe], label = r"Along "+axe)
    for axe in ax:
        axe.legend()
    return fig





    
    

    


    

