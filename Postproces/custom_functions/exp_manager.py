# from Info_making import read_info
from IPython.display import display, Markdown
from custom_functions.read import read_json, read_observables
import glob
import scipy.constants as cst
from custom_functions.fouriertransform import *
from custom_functions.high_harmonics import findHarmPeaks
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

class Experience:
    '''A class for and EDUS experience.

    Parameters
    ----------
    path : str
        Path to the results folder returned by EDUS. Folder must have the .txt files and the .json input file OR the .json input file and a subfolder /Outputs containing the .txt files.
    laser_index : int, optional
        Index of the laser to consider as main. The HHG fundamental frequency will be the same as this laser.
    
    Attributes
    ----------
    path : str, private
        Path to the results folder returned by EDUS. Folder must have the .txt files and the .json input file OR the .json input file and a subfolder /Outputs containing the .txt files.
    number : int, global
        ID of the instance.

    tb_file : str
        Name of the tb file used to run the EDUS calculation.
    filled_bands : int
        Number of filled bands.
    dt : float 
        Time step of the propagation of EOMs by EDUS.
    total_time : float
        Total simulation time, in fs.
    grid : list of int
        K-points grid used in the EDUS calculation.
    lasers : list of Lasers
        Lasers used in the EDUS calculation.
    main_laser : 
        Instance of Laser considered to be the driving one.
    
    time_au : array
        Array of all the time steps in atomic units.
    population_time : array
        Array of the population evolution in temporal domain.
    field_time_au : array
        Array of the applied field, in atomic units, in time domain.
    velocity_time_au : array
        Array of the velocity, in atomic units, in time domain. 
    field_freq_au : array
        Array of the applied field, in atomic units, in frequency domain.
    frequ_axis_eV : array
        Array of the energy in the frequency domain, in electronvolts.
    frequ_axis_au : array
        Array of the energy in the frequency domain, in atomic units.

    HHG_fund_freq_eV : float
        Fundamental frequency, in eV, of the HHG.
    HHG_acceleration_freq_au : array
        Norm of the acceleration, in atomic units and in the frequency domain.

    HHG_harm_idx : array
        Indices of the maximas of the HHG spectra.
    HHG_pics_height : array
        Value of the maximas of the HHG sepctra.
    HHG_normed_pics_height : array
        Value of the maximas of the HHG sepctra, normalized by the first pic intensity.

    velocity_FT_dict : dict of arrays
        Stored arrays of the fourier transform (FT) of the velocity, using a defined smearing. Keys are :code:`f'{smearing}'`.
    HHG_spectrum : dict of HHG_Spectra
        Dictionary storing all the HHG spectras computed, using the parameters from :code:`gen_HHG_spectra(**kwargs)`method. Keys are in the form :code:`f"{key}_{value}" for key, value in **kwargs.items`.
    HHG_peaks_spectrum : dict of HHG_Peaks_Spectra
        Dictionary storing all the HHG peaks spectras computed, using the parameters from :code:`gen_HHG_peaks_spectra(**kwargs)`method. Keys are in the form :code:`f"{key}_{value}" for key, value in **kwargs.items`.

    Methods
    ---------
    `print_info()`
        For jupyter notbooks usage, will pring a md recap of useful information of the experience.
    
    `gen_observables(smearing=0.)`
        If not already done, generates the observables of the experience, as well as the FT of the field and the velocity using the provided smearing (only for velocity)
    `gen_HHG_data(smearing=0.)`
        If not already done, computes the FT of the acceleration from the FT of the velocity for the given smearing.
    `gen_HHG_peaks(max_harmonic=20, smearing=0.)`
        Finds the height of the emission at every harmonic order from 1 up to max_harmonics. The emission spectra is computed from :code:`gen_HHG_data(smearing)`.
    `gen_HHG_spectra(start=0,stop=10,axes=[],cut_eV=0.)`
        Generates an instance of :code:`HHG_Spectra(**kwargs)` and stores it in :code:`HHG_spectrum`.
    `gen_HHG_peaks_spectra(max_harmonic=20,normed = False, axes = [])`
        Generates an instance of :code:`HHG_Peak_Spectra(**kwargs)` and stores it in :code:`HHG_peaks_spectrum`.

    `del_observables()`
        Deletes all stored observables previously generated with :code:`gen_observables`.
    `del_figures()`
        Deletes all stored figures.
    `del_HHG_peaks()`
        Deletes peaks data. 

    `listFigures()`
        List all the figures stored.
    `getFigure(key)`
        Returns the :code:`HHG_Spectra` or :code:`HHG_Peaks_Spectra` instance corresponding to the provided key. 
    '''

    Exp_number = 0

    @property
    def path(self):
        return self.__path
    
    def __init__(self, path, laser_index =0):
        Experience.Exp_number += 1
        self.number = Experience.Exp_number
        self.__path = path 
        self.__laser_index = laser_index
        info_dict = read_info(path)

        self.tb_file = info_dict['tb_file'].split('/')[-1] + '_tb.dat'
        self.filled_bands = info_dict['filledbands']
        lasers_list = info_dict['lasers']
        self.dt = info_dict['dt']
        self.total_time = info_dict['finaltime']
        self.grid = info_dict['grid']
        self.lasers = []
        for laser in lasers_list:
            self.lasers.append(Laser(self,laser['intensity'], laser['frequency'], laser['frequency_units'],laser['polarization'], laser['cycles'], laser['t0']))
        
        self.velocity_FT_dict = {}

        self.HHG_spectrum = {}
        self.HHG_peaks_spectrum = {}

        self.main_laser = self.lasers[self.__laser_index]

        self.gen_observables()

        ###### Generate the info string ###### 
        self._info = f"## Experience n° {self.number}\n"

        formatted_dt = f"${self.dt / 10**int(f'{self.dt:e}'.split('e')[1]):.2f} \\cdot" + r"10^{" + f"{{{int(f'{self.dt:e}'.split('e')[1])}}}" + r"}$"

        self._info = self._info + f'### Simulation info:\n* Hamiltonian file : {self.tb_file} \n* Number of filled bands : {self.filled_bands}\n* Time step : {formatted_dt} fs\n* Total time : ${self.total_time}$ fs\n* Grid : ${self.grid[0]}\\times {self.grid[1]}$'

        for laser in self.lasers:

            formatted_i = f"${laser.intensity / 10**int(f'{laser.intensity:e}'.split('e')[1]):.2f} \\cdot" + r"10^{" + f"{{{int(f'{laser.intensity:e}'.split('e')[1])}}}" + r"}$"

            formatted_freq = f"${laser.frequ_Hz / 10**int(f'{laser.frequ_Hz:e}'.split('e')[1]):.3f} \\cdot" + r"10^{" + f"{{{int(f'{laser.frequ_Hz:e}'.split('e')[1])}}}" + r"}$"

            self._info = self._info + '\n ### Laser info:' f'\n* Frequency : ${laser.frequ_eV}$ eV *i.e.* {formatted_freq} Hz' + f'\n * Intensity : {formatted_i}'+  ' w.cm $^{-2}$' + f'\n * Cycles : ${laser.cycles}$' + f'\n * Total pulse duration : ${laser.total_pulse_time:.2f}$ fs'

        
        

    def __str__(self):
        return f"Experience n° {self.number}"
        

    def gen_observables(self, smearing = 0.):
        if not hasattr(self, 'time_au'):
            try:
                self.time_au, self.population_time, self.field_time_au, self.velocity_time_au = read_observables(self.path, 'new')
            except FileNotFoundError:
                self.time_au, self.population_time, self.field_time_au, self.velocity_time_au = read_observables(self.path + '/Output', 'new')
            self.frequ_axis_eV, self.field_freq_au = FourierTransform(self.time_au, self.field_time_au)
            self.frequ_axis_au = self.frequ_axis_eV*constants.physical_constants["hartree-electron volt relationship"][0]

        if f'{smearing}' in self.velocity_FT_dict.keys():
            pass
        else:
            _, velocity_frequ_au = FourierTransform(self.time_au, self.velocity_time_au, cut_eV = smearing)
            
            self.velocity_FT_dict[f'{smearing}'] = velocity_frequ_au
            
    def gen_HHG_data(self, smearing = 0.):
        if f'{smearing}' not in self.velocity_FT_dict.keys():
            self.gen_observables(smearing=smearing)

        self.HHG_fund_freq_eV = self.lasers[self.__laser_index].frequ_eV

        Aw_au = self.velocity_FT_dict[f'{smearing}'] * 1j * self.frequ_axis_au

        self.HHG_acceleration_freq_au_x = Aw_au[0,:]
        self.HHG_acceleration_freq_au_y = Aw_au[1,:]
        self.HHG_acceleration_freq_au_z = Aw_au[2,:]

        self.HHG_acceleration_freq_au = np.linalg.norm(Aw_au, axis=0) #norm over all axis


    def gen_HHG_peaks(self,max_harmonic=20, smearing = 0.): ##### OBSOLETE SINCE NEW PEAKS SEARCH ! 

        self.gen_HHG_data(smearing=smearing) #generating HHG data is fast, can just select the one needed.

        self.HHG_harm_idx = findHarmPeaks(self.frequ_axis_eV/self.HHG_fund_freq_eV, self.HHG_acceleration_freq_au, limit=max_harmonic)
        self.HHG_pics_height = self.HHG_acceleration_freq_au[self.HHG_harm_idx]**2
        self.HHG_normed_pics_height = self.HHG_pics_height**2/np.max(self.HHG_pics_height)**2

        ######## No Axes Implementation for now ###########

        # temp_dict['HHG_harm_idx_x'] = findHarmPeaks(self.frequ_axis_eV/self.HHG_fund_freq_eV, np.abs(self.HHG_acceleration_freq_au_x), limit=max_harmonic, cutoff=cutoff)
        # temp_dict['HHG_harm_idx_y'] = findHarmPeaks(self.frequ_axis_eV/self.HHG_fund_freq_eV, np.abs(self.HHG_acceleration_freq_au_y), limit=max_harmonic, cutoff=cutoff)
        # temp_dict['HHG_harm_idx_z'] = findHarmPeaks(self.frequ_axis_eV/self.HHG_fund_freq_eV, np.abs(self.HHG_acceleration_freq_au_z), limit=max_harmonic, cutoff=cutoff)



    def gen_HHG_spectra(self, start = 0, stop=10, axes=[], smearing = 0.):
        if f'{smearing}' not in self.velocity_FT_dict.keys():
            print("Generating observables.")
            self.gen_observables(smearing=smearing)
        
        self.gen_HHG_data()
        # if f'Start_{start}-Stop_{stop}-Axes_{axes}-Cut_eV_{cut_eV}' in self.HHG_spectrum.keys():
        #     print(f"Spectra already computed, plotting it now with Experience_{self.number}.plot_HHG_spectra(*args)\nYou can access it with Experience_{self.number}.HHG_spectrum['Start_{start}-Stop_{stop}-Axes_{axes}-Cut_eV_{cut_eV}']")
        #     self.plot_HHG_spectra(start=start, stop=stop, axes=axes, cut_eV=cut_eV)
        # else:
            # self.HHG_spectrum[f'Start_{start}-Stop_{stop}-Axes_{axes}-Cut_eV_{cut_eV}'] = HHG_Spectra(self, start=start, stop= stop, axes=axes, cut_eV=cut_eV)
        self.HHG_spectrum[f'Start_{start}-Stop_{stop}-Axes_{axes}-smearing_{smearing}'] = HHG_Spectra(self, start=start, stop= stop, axes=axes, smearing=smearing)

    def gen_HHG_peaks_spectra(self, max_harmonic=20,normed = False, axes = []):
        if not hasattr(self, 'time_au'):
            print("Generating observables.")
            self.gen_observables()
        if not hasattr(self, 'HHG_fund_freq_eV'):
            print("Generating HHG data.")
            self.gen_HHG_data()
        if not hasattr(self, 'HHG_harm_idx'):
            print("Generating peaks indexes.")
            self.gen_HHG_peaks(max_harmonic=max_harmonic)
        
        if f'Max_harm_{max_harmonic}-Normed_{normed}-Axes_{axes}' in self.HHG_peaks_spectrum.keys():
            print(f"Spectra already computed, plotting it now with Experience_{self.number}.plot_HHG_peaks_spectra(*args)\nYou can access it with Experience_{self.number}.HHG_peaks_spectrum['Max_harm_{max_harmonic}-Normed_{normed}-Axes_{axes}']")
            self.plot_HHG_peaks_spectra(max_harmonic=max_harmonic, axes=axes, normed=normed)
        else:
            self.HHG_peaks_spectrum[f'Max_harm_{max_harmonic}-Normed_{normed}-Axes_{axes}'] = HHG_Peak_Spectra(self, max_harmonic=max_harmonic, axes=axes, normed=normed)
        # self.HHG_peaks_spectrum[f'Cutoff_{cutoff}-Max_harm_{max_harmonic}-Normed_{normed}-Axes_{axes}'] = HHG_Peak_Spectra(self, max_harmonic=max_harmonic, axes=axes, cutoff=cutoff, normed=normed)
        
    def del_figures(self):
        self.HHG_spectrum = {}
        self.HHG_peaks_spectrum = {}

    def del_observables(self):
        self.velocity_FT_dict = {}
        self.time_au, self.population_time, self.field_time_au, self.velocity_time_au = None, None, None,None
        self.frequ_axis_eV, self.field_freq_au = None, None
        self.frequ_axis_au = None

    def del_HHG_peaks(self):
        self.HHG_harm_idx = None
        self.HHG_pics_height = None

    
    def print_info(self):
        display(Markdown(self._info))

    def plot_HHG_spectra(self, start = 0, stop=10, axes=[], smearing = 0.):
        """Plots the HHG spectra with the specified parameters

        Returns
        -------
        Spectra : HHG_Spectra
        """
        if f'Start_{start}-Stop_{stop}-Axes_{axes}-smearing_{smearing}' not in self.HHG_spectrum.keys():
            self.gen_HHG_spectra(start = start, stop=stop, axes=axes, smearing = smearing)
            return self.HHG_spectrum[f'Start_{start}-Stop_{stop}-Axes_{axes}-smearing_{smearing}']
        else:
            self.HHG_spectrum[f'Start_{start}-Stop_{stop}-Axes_{axes}-smearing_{smearing}'].show()
            return self.HHG_spectrum[f'Start_{start}-Stop_{stop}-Axes_{axes}-smearing_{smearing}']

    def plot_HHG_peaks_spectra(self, max_harmonic=20,normed = False, axes = []):
        if f'Max_harm_{max_harmonic}-Normed_{normed}-Axes_{axes}' not in self.HHG_peaks_spectrum.keys():
            self.gen_HHG_peaks_spectra( max_harmonic=max_harmonic,normed = normed, axes = axes)
            return self.HHG_peaks_spectrum[f'Max_harm_{max_harmonic}-Normed_{normed}-Axes_{axes}']
        else:
            self.HHG_peaks_spectrum[f'Max_harm_{max_harmonic}-Normed_{normed}-Axes_{axes}'].show()
            return self.HHG_peaks_spectrum[f'Max_harm_{max_harmonic}-Normed_{normed}-Axes_{axes}']

    def listFigures(self):
        fig_list = []
        figures = "Available figures:\n  HHG Spectras:\n"

        for figure in self.HHG_spectrum.keys():
            figures = figures + f'       * {figure}\n'
            fig_list.append(figure)
        
        figures = figures + '   HHG Peaks Spectras:\n'

        for figure in self.HHG_peaks_spectrum.keys():
            figures = figures + f'       * {figure}\n'
            fig_list.append(figure)

        print(figures)
        return fig_list

    def getFigure(self,key):
        if key in self.HHG_spectrum.keys():
            return self.HHG_spectrum[key]
        elif key in self.HHG_peaks_spectrum.keys():
            return self.HHG_peaks_spectrum[key]
        else:
            raise KeyError("Figure not found, please check listFigures()")





class Laser:
    """A class for the lasers used in a EDUS simulation. Note that all the parameters are also attributes, except :code:`experience` and :code:`frequency`.

    Parameters
    ----------
    experience : instance of Experience
        The parent experience of the laser.
    frequency : float
        Frequency of the laser.
    intensity : float
        Intensity of the laser, in w/cm^2
    frequ_units : str
        Unit of the laser's frequency.
    polarization : list of floats
        Polarization of the laser in the x,y,z directions.
    cycles : int
        Number of optical cycles.
    t0 : float
        Starting time, in femtoseconds, of the laser pulse.

    Attributes
    ----------
    frequ_Hz : float
        Laser's frequency in hertz.
    frequ_eV : float
        Laser's frequency in electronvolts.
    total_pulse_time : float
        Total time lenght of the laser.

    Methods
    -------
    `plot(time_window = (0,0), axes=['x','y','z'])`
        Plots the laser pulses and returns the associated figure.
    
    """
    def __init__(self,experience, intensity, frequency, frequ_units, polarization, cycles, t0):
        self.intensity = intensity
        self.frequ_units = frequ_units
        self.polarization = polarization
        self.cycles = cycles
        self.t0 = t0
        self.__experience = experience

        if frequ_units == 'electronvolt':
            self.frequ_eV = frequency
        else:
            pass # WIP

        self.frequ_Hz = self.frequ_eV*cst.physical_constants['electron volt-hertz relationship'][0]

        self.total_pulse_time = (1/self.frequ_Hz) * 1e15 * self.cycles

    def plot(self, time_window = (0,0), axes=['x','y','z']):
        """Plots the laser pulses.

    Parameters
    ----------
    datadir : str
        Path to data folder.
    time_window : tuple of floats, optional
        Time interval of the simulation to consider in au. If none provided, the whole laser is plotted.
    axes : list of str, optional
        Direction of polarisation of the laser to plot. If left empty, will plot the whole norm of the acceleration.

    Returns
    -------
    fig : matplotlib.pyplot.Figure
        The absorption spectra.

    Raises
    -------
    ValueError
        If the time_window provided is not in proper order.
    ValueError
        If the axes are not provided properly, i.e. different than x, y and z or empty.
    """
        
        if time_window[0] > time_window[1]:
            raise ValueError(f"You cannot reverse time !")
    
        if ('x' not in axes) and ('y' not in axes) and ('z' not in axes):
            raise ValueError(f"Axis are poorly defined")
        if not hasattr(self.__experience, 'time_au'):
            self.__experience.gen_observables()

        time_au = self.__experience.time_au
        laser_au = self.__experience.field_time_au

         ###### Select proper time window ######
        if time_window != (0,0):
            idx = np.where(np.logical_and(time_au > time_window[0], time_au < time_window[1]))[0]
            time_au = time_au[idx[0]:idx[-1]]
            laser_au = laser_au[:,idx[0]:idx[-1]]
            print("First time step: ", time_au[0],"a.u. \nLast time step: ", time_au[-1], "a.u.")

        laserA_au_x= laser_au[0,:] 
        laserA_au_y =  laser_au[1,:] 
        laserA_au_z = laser_au[2,:]

        ###### Plotting ######
        fig, ax = plt.subplots()
        
        for axe in axes:
            ax.plot(time_au, locals()["laserA_au_" + axe], label = r"Along " + axe)
        
        ax.legend()
        return fig
        



class HHG_Spectra:
    '''A class for high harmonic generation spectra generated from an Experience. 

    Parameters
    ----------
    Experience : instance of Experience
        The considered experience, the observables must be previouly computed with Experience.gen_HHG_data().
    start : int, optional
        Starting harmonic order to plot.
    stop : int, optional
        Highest harmonic order to plot.
    axes : list of str, optional
        Axes to plot.
    '''
    def __init__(self, Experience, start = 0, stop=10, axes=[], smearing = 0.):

        if ('x' not in axes) and ('y' not in axes) and ('z' not in axes) and axes!=[]:
            raise ValueError(f"Axis are poorly defined.")
    

        self.figure, self.ax = plt.subplots()

        self.exp_number = Experience.number
        self.start = start
        self.stop = stop
        self.axes = axes
        self.smearing = smearing


        x_harmonic_axis = Experience.frequ_axis_eV/Experience.HHG_fund_freq_eV

        if axes == []:
            self.ax.plot(x_harmonic_axis, Experience.HHG_acceleration_freq_au**2)

        for axe in axes:
            self.ax.plot(x_harmonic_axis, np.abs(vars(Experience)["HHG_acceleration_freq_au_"] + axe)**2, label = r"Along " + axe)

        self.ax.set_yscale('log')

        self.ax.set_xticks(np.linspace(start,stop, stop +1-start, dtype = int)) # get all the orders
        self.ax.set_xlim(start,stop)
        self.ax.set_xlabel(r"Harmonic order")
        self.ax.set_ylabel(r"Intensity (arb. units, log scale)")
        self.ax.legend()

        self.__base_ylims = self.ax.set_ylim()


    def __str__(self):
        form_ax = formatAxes(self.axes)
        return f'HHG spectra from experience n° {self.exp_number} with parameters:\nFirst harmonic: {self.start}\nLast harmonic: {self.stop}\nAxes: {form_ax}'

    def show(self):
        plotFigures(self.figure)
    def get_axes(self):
        return self.ax
    
    def resetAxes(self):
        self.ax.set_yscale('log')
        self.ax.set_ylim(self.__base_ylims[0],self.__base_ylims[1])
        self.ax.set_xticks(np.linspace(self.start,self.stop, self.stop +1-self.start, dtype = int))
        self.ax.set_xlim(self.start,self.stop)
        self.ax.set_xlabel(r"Harmonic order")
        self.ax.set_ylabel(r"Intensity (arb. units, log scale)")
        self.ax.set_title(r"")
        self.ax.grid(False)
        self.ax.legend()


class HHG_Peak_Spectra():
    def __init__(self, Experience, max_harmonic=20, axes=[], cutoff = 1.5e-1, normed = False):
        self.figure, self.ax = plt.subplots()

        self.exp_number = Experience.number
        self.max_harmonic = max_harmonic
        self.axes = axes
        self.cutoff = cutoff
        self.normed = normed

        x_harmonic_axis = Experience.frequ_axis_eV/Experience.HHG_fund_freq_eV 

        harm_idx = Experience.HHG_harm_idx

        # harm_idx_x = Experience.HHG_peaks_dict[f'{cutoff}']['HHG_harm_idx_x']
        # harm_idx_y = Experience.HHG_peaks_dict[f'{cutoff}']['HHG_harm_idx_y']   
        # harm_idx_z = Experience.HHG_peaks_dict[f'{cutoff}']['HHG_harm_idx_z']   



        if axes == []:
            if normed:
                y_axis = Experience.HHG_normed_pics_height
            else:
                y_axis = Experience.HHG_pics_height
            self.ax.scatter(x_harmonic_axis[harm_idx], y_axis)

        for axe in axes:
            if normed:
                y_axis = np.abs(vars(Experience)["HHG_acceleration_freq_au_"+ axe] [vars()['harm_idx'+axe]])**2/np.max(np.abs(vars(Experience)["HHG_acceleration_freq_au_"+ axe] [vars()['HHG_harm_idx'+axe]])**2)
            else:
                y_axis = np.abs(vars(Experience)["HHG_acceleration_freq_au_"+ axe] [vars()['harm_idx'+axe]])**2
            self.ax.plot(x_harmonic_axis[vars()['harm_idx'+axe]] , y_axis, label = r"Along " + axe)

        self.ax.set_yscale('log')

        self.ax.set_xticks(np.linspace(0,max_harmonic, max_harmonic +1, dtype = int)) # get all the orders
        self.ax.set_xlim(0,max_harmonic)
        self.ax.set_xlabel(r"Harmonic order")
        self.ax.set_ylabel(r"Intensity (arb. units, log scale)")
        self.ax.legend()

        self.__base_ylims = self.ax.set_ylim()

    
    def __str__(self):
        form_ax = formatAxes(self.axes)
        return f'HHG pics spectra from experience n° {self.exp_number} with parameters:\nMax harmonic: {self.max_harmonic}\nWindow to detect peak: {self.cutoff}\nAxes: {form_ax}\nNormed : f{self.normed}'

    def show(self):
        plotFigures(self.figure)

    def get_axes(self):
        return self.ax
    def resetAxes(self):
        self.ax.set_yscale('log')
        self.ax.set_xticks(np.linspace(0,self.max_harmonic, self.max_harmonic +1, dtype = int))
        self.ax.set_xlim(0,self.max_harmonic)
        self.ax.set_ylim(self.__base_ylims[0],self.__base_ylims[1])
        self.ax.set_xlabel(r"Harmonic order")
        self.ax.set_ylabel(r"Intensity (arb. units, log scale)")
        self.ax.set_title(r"")
        self.ax.grid(False)
        self.ax.legend()

def read_info(datadir):
    return read_json(glob.glob(datadir + '/*.json')[0])

def formatAxes(axes_list):
    if axes_list == []:
        return "Normed"
    formatted = ''
    if 'x' in axes_list:
        formatted = formatted + ' x,'
    if 'y' in axes_list:
        formatted = formatted + ' y,'
    if 'z' in axes_list:
        formatted = formatted + ' z,'
    return formatted


def plotFigures(fig):
    axis = fig.get_axes()
    if len(np.array(axis).shape) == 0:
            figure, axes = plt.subplots()
            for line in axis.get_lines():
                axes.plot(line.get_xdata(), line.get_ydata(), label = line.get_label(), color= line.get_color(), linestyle = line.get_linestyle())
            for collection in axis.collections:
                offsets = collection.get_offsets()
                axes.scatter(offsets[:,0], offsets[:,1])
            for text in axis.texts:  # Handle text elements
                axes.text(text.get_position()[0], text.get_position()[1], text.get_text(), fontsize=text.get_fontsize(), color=text.get_color())
                
            axes.set_xlim(ax.get_xlim())
            axes.set_ylim(ax.get_ylim())

            axes.set_xscale(ax.get_xscale())
            axes.set_yscale(ax.get_yscale())

            axes.set_xticks(ax.get_xticks())
            axes.set_yticks(ax.get_yticks())

            axes.set_xlabel(ax.get_xlabel())
            axes.set_ylabel(ax.get_ylabel())
            axes.set_title(ax.get_title())

            axes.grid(ax.get_xgridlines()[0].get_visible())

            axes.legend()


    if len(np.array(axis).shape) == 1:
        try:
            size = np.array(axis).shape[1]
        except IndexError:
            size = np.array(axis).shape[0]

        figure, axes = plt.subplots(size)
        if size == 1:
            axes = [axes]
        for i in range(size):
            ax = axis[i]
            for line in ax.get_lines():
                axes[i].plot(line.get_xdata(), line.get_ydata(), label = line.get_label(), color= line.get_color(), linestyle = line.get_linestyle())
            for collection in ax.collections:
                offsets = collection.get_offsets()
                axes[i].scatter(offsets[:,0], offsets[:,1])
            for text in ax.texts:  # Handle text elements
                axes[i].text(text.get_position()[0], text.get_position()[1], text.get_text(), fontsize=text.get_fontsize(), color=text.get_color())

            axes[i].set_xlim(ax.get_xlim())
            axes[i].set_ylim(ax.get_ylim())

            axes[i].set_xscale(ax.get_xscale())
            axes[i].set_yscale(ax.get_yscale())

            axes[i].set_xticks(ax.get_xticks())
            axes[i].set_yticks(ax.get_yticks())

            axes[i].set_xlabel(ax.get_xlabel())
            axes[i].set_ylabel(ax.get_ylabel())
            axes[i].set_title(ax.get_title())

            axes[i].grid(ax.get_xgridlines()[0].get_visible())
            axes[i].legend()

    if len(np.array(ax).shape) == 2:
        figure, axes = plt.subplots(np.array(axis).shape[0],np.array(axis).shape[1])

        for i, axI in enumerate(axis):
            for j, ax in enumerate(axI):
                for line in ax.get_lines():
                    axes[i][j].plot(line.get_xdata(), line.get_ydata(), label = line.get_label(), color= line.get_color(), linestyle = line.get_linestyle())
                for collection in ax.collections:
                    offsets = collection.get_offsets()
                    axes[i][j].scatter(offsets[:,0], offsets[:,1])
                for text in ax.texts:  # Handle text elements
                    axes[i][j].text(text.get_position()[0], text.get_position()[1], text.get_text(), fontsize=text.get_fontsize(), color=text.get_color())

                axes[i][j].set_xlim(ax.get_xlim())
                axes[i][j].set_ylim(ax.get_ylim())

                axes[i][j].set_xscale(ax.get_xscale())
                axes[i][j].set_yscale(ax.get_yscale())

                axes[i][j].set_xticks(ax.get_xticks())
                axes[i][j].set_yticks(ax.get_yticks())

                axes[i][j].set_xlabel(ax.get_xlabel())
                axes[i][j].set_ylabel(ax.get_ylabel())
                axes[i][j].set_title(ax.get_title())

                axes[i][j].grid(ax.get_xgridlines()[0].get_visible())

                axes[i][j].legend()
    plt.show()
    return figure

def joinFigures(fig1, fig2, samefig = False):
    """Joins two matplotlib.pyplot.Figure into a new one.

    Parameters
    ----------
    fig1 : matplotlib.pyplot.Figure
        First figure to join.
    fig2 : matplotlib.pyplot.Figure
        Second figure to join.

    Returns
    -------
    mergedFigure : matplotlib.pyplot.Figure
        A figure containing all the axes of the figures. 
    """

    fig1Axes = fig1.get_axes()
    fig2Axes = fig2.get_axes()

    try:
        len(fig1Axes)
    except TypeError:
        fig1Axes = [fig1Axes]
    try:
        len(fig2Axes)
    except TypeError:
        fig2Axes = [fig2Axes]
    

    if samefig:
        mergedFigure, axes = plt.subplots()
        axes = [axes]
    else:
        mergedFigure, axes = plt.subplots(len(fig1Axes) + len(fig2Axes), 1)

    i = 0 
    for ax in fig1Axes:
        for line in ax.get_lines():
            axes[i].plot(line.get_xdata(), line.get_ydata(), label = line.get_label(), color= line.get_color(), linestyle = line.get_linestyle())
        for collection in ax.collections:
            offsets = collection.get_offsets()
            axes[i].scatter(offsets[:,0], offsets[:,1])
        for text in ax.texts:  # Handle text elements
            axes[i].text(text.get_position()[0], text.get_position()[1], text.get_text(), fontsize=text.get_fontsize(), color=text.get_color())

        axes[i].set_xlim(ax.get_xlim())
        axes[i].set_ylim(ax.get_ylim())

        axes[i].set_xscale(ax.get_xscale())
        axes[i].set_yscale(ax.get_yscale())

        axes[i].set_xticks(ax.get_xticks())
        axes[i].set_yticks(ax.get_yticks())

        axes[i].set_xlabel(ax.get_xlabel())
        axes[i].set_ylabel(ax.get_ylabel())
        axes[i].set_title(ax.get_title())
        axes[i].legend()
        if not samefig:
            i+=1
    for ax in fig2Axes:
        for line in ax.get_lines():
            axes[i].plot(line.get_xdata(), line.get_ydata(), label = line.get_label(), color= line.get_color(), linestyle = line.get_linestyle())
        for collection in ax.collections:
            offsets = collection.get_offsets()
            axes[i].scatter(offsets[:,0], offsets[:,1])
        for text in ax.texts:  # Handle text elements
            axes[i].text(text.get_position()[0], text.get_position()[1], text.get_text(), fontsize=text.get_fontsize(), color=text.get_color())

        axes[i].set_xlim(ax.get_xlim())
        axes[i].set_ylim(ax.get_ylim())

        axes[i].set_xscale(ax.get_xscale())
        axes[i].set_yscale(ax.get_yscale())

        axes[i].set_xticks(ax.get_xticks())
        axes[i].set_yticks(ax.get_yticks())

        axes[i].set_xlabel(ax.get_xlabel())
        axes[i].set_ylabel(ax.get_ylabel())
        axes[i].set_title(ax.get_title())
        axes[i].legend()
        if not samefig:
            i+=1

    return mergedFigure
