import json
import numpy as np
from custom_functions.read import read_json
import glob
import scipy.constants as cst
from IPython.display import display, Markdown

def make_info(datadir):
    data = read_info(datadir)
    info = ""

    tb_file = data['tb_file'].split('/')[-1] + '_tb.dat'
    filled_bands = data['filledbands']

    lasers = data['lasers']
    dt = data['dt']

    formatted_dt = f"${dt / 10**int(f'{dt:e}'.split('e')[1]):.2f} \\cdot" + r"10^{" + f"{{{int(f'{dt:e}'.split('e')[1])}}}" + r"}$"

    grid = data['grid']
    total_time = data['finaltime']
    info = info + f'### Simulation info:\n* Hamiltonian file : {tb_file} \n* Number of filled bands : {filled_bands}\n* Time step : {formatted_dt} fs\n* Total time : ${total_time}$ fs\n* Grid : ${grid[0]}\\times {grid[1]}$'

    for laser in lasers:
        intensity = laser['intensity']
        frequency = laser['frequency']
        cycles = laser['cycles']

        hz_frequ = frequency*cst.physical_constants['electron volt-hertz relationship'][0]
        
        total_pulse_time = (1/hz_frequ) * 1e15 * cycles
        
        
        formatted_i = f"${intensity / 10**int(f'{intensity:e}'.split('e')[1]):.2f} \\cdot" + r"10^{" + f"{{{int(f'{intensity:e}'.split('e')[1])}}}" + r"}$"

        formatted_freq = f"${hz_frequ / 10**int(f'{hz_frequ:e}'.split('e')[1]):.3f} \\cdot" + r"10^{" + f"{{{int(f'{hz_frequ:e}'.split('e')[1])}}}" + r"}$"


        info = info + '\n ### Laser info:' f'\n* Frequency : ${frequency}$ eV *i.e.* {formatted_freq} Hz' + f'\n * Intensity : {formatted_i}'+  ' w.cm $^{-2}$' + f'\n * Cycles : ${cycles}$' + f'\n * Total pulse duration : ${total_pulse_time:.2f}$ fs'


    display(Markdown(info))


def read_info(datadir):
    return read_json(glob.glob(datadir + '/*.json')[0])

def gen_info(info_dict, key_list):

    for key in key_list:
        if key not in info_dict.keys():
            raise ValueError(f"The provided key is not valid: {key}")
        
    values = []

    for key in key_list:
        vars()[key] = info_dict[key]
    

    # dict_ = {}

    # dict_["tb_file"]           = js["tb_file"]
    # dict_["coulomb"]           = js["coulomb"]
    # dict_["grid"]              = js['grid']
    # dict_["filledbands"]       = js["filledbands"]
    # dict_["dt"]                = js["dt"]
    # dict_["printresolution"]   = js["printresolution"]
    # dict_["initialtime"]       = js["initialtime"]
    # dict_["finaltime"]         = js["finaltime"]
    # dict_["num_bands"]         = js["num_bands"]

    # dict_["lasers"]            = 


    # tb_file = data['tb_file'].split('/')[-1] + '_tb.dat'
    # filled_bands = data['filledbands']
    # lasers = data['lasers']
    # dt = data['dt']
    