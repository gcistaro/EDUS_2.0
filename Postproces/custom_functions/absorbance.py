
import numpy as np
from custom_functions.fouriertransform import FourierTransform
from scipy import constants
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.io as pio

def absorbance(t_au, Vt_au, Et_au, limits=[]):
    print("Fourier transforming")
    freq_eV, Jw_au = FourierTransform(t_au, Vt_au, True)
    _, Ew_au       = FourierTransform(t_au, Et_au, False)
    print("done")


    #cut a window
    if limits == []:
        Ew_au_abs_max = np.max(np.abs(Ew_au))
        Ew_au_condition = np.zeros(len(freq_eV))
        for i in range(3):
            Ew_au_condition += np.abs(np.abs(Ew_au[i,:]) >= 0.05*Ew_au_abs_max)
        freq_eV_condition = freq_eV >= 0.0
        meaningful_Ew_au_condition = freq_eV_condition*Ew_au_condition
        index = np.where(meaningful_Ew_au_condition != 0)[0]
    else:
        index = np.where(np.logical_and( freq_eV > limits[0], freq_eV < limits[1] ))[0]

    plt.plot(freq_eV, np.abs(Ew_au[0, :]))
    plt.plot(freq_eV[index], np.abs(Ew_au[0, index]))
    plt.xlim(left=0)
    plt.show()
    plt.close()

    Ew_au_fig = go.Figure()
    Ew_au_fig.add_trace(go.Scatter(x=freq_eV, y=np.abs(Ew_au[0, :]), mode='lines'))
    Ew_au_fig.add_trace(go.Scatter(x=freq_eV[index], y=np.abs(Ew_au[0, index]), mode='lines'))
    Ew_au_fig.update_layout(xaxis=dict(range=[0, 2*max(freq_eV[index])]))
    pio.write_html(Ew_au_fig, 'Ew_au.html')

    Jw_au = Jw_au[:,index]
    Ew_au = Ew_au[:,index]
    freq_eV = freq_eV[index]
    freq_au = freq_eV*constants.physical_constants["hartree-electron volt relationship"][0]

    freq_au = freq_eV*constants.physical_constants["hartree-electron volt relationship"][0]
    alpha = constants.physical_constants["fine-structure constant"][0]

    linear_conductivity = Jw_au[0]/Ew_au[0]
    np.savetxt('linear_conductivity.txt', np.transpose([freq_eV, linear_conductivity]))
    plt.figure(figsize=(10, 7.5))
    plt.plot(freq_eV, np.real(linear_conductivity), label = r'Real{$\omega$}')
    plt.plot(freq_eV, np.imag(linear_conductivity), label = r'Imag{$\omega$}')
    plt.xlabel(r'$\omega$ (eV)')
    plt.ylabel(r'$\tilde{\sigma}}(\omega)$')
    plt.legend()
    plt.savefig('linear_conductivity', dpi=800)
    plt.show()
    plt.close()

    dw_au = -1j*Jw_au/freq_au # constants are missing (volume)
    Sw_au  = (2*np.imag( np.sum( [ np.conj(Ew_au[ix])*dw_au[ix] for ix in [0,1,2] ], axis=0) ) )        #response function
    Iw_au = 2*np.sum( [np.abs(Ew_au[ix])*np.abs(Ew_au[ix]) for ix in [0,1,2] ], axis=0)/(4.*np.pi*alpha*freq_au) # not sure about these constants
    plt.plot(freq_eV, np.real(dw_au[0]), label = r'Real($\tilde{\left\langle \mathbf{d}\right\rangle }(\omega)$)')
    plt.plot(freq_eV, np.imag(dw_au[1]), label = r'Imag($\tilde{\left\langle \mathbf{d}\right\rangle }(\omega)$)')
    plt.legend()
    plt.show()
    
    print("Ew[1] min and max: ", np.min(np.abs(Ew_au[1])), np.max(np.abs(Ew_au[1])))
    print("Ew[2] min and max: ", np.min(np.abs(Ew_au[2])), np.max(np.abs(Ew_au[2])))
    print("Ew min and max: ", np.min(np.abs(Ew_au)), np.max(np.abs(Ew_au)))
    print("Vw min and max: ", np.min(np.abs(Jw_au)), np.max(np.abs(Jw_au)))
    print("rw min and max: ", np.min(np.abs(dw_au)), np.max(np.abs(dw_au)))
    print("Sw min and max: ", np.min(np.abs(Sw_au)), np.max(np.abs(Sw_au)))
    print("Iw min and max: ", np.min(np.abs(Iw_au)), np.max(np.abs(Iw_au)))
    Absorbance =Sw_au/Iw_au

    return freq_eV, Absorbance
