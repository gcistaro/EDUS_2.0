import numpy as np
import os
from scipy import constants


def read_observables(folder, version, print_info = False):
    ''' Reads folder and extracts Time, Population, Laser and Velocity

    Returns
    -------
    t_au
        Time 
    pt
        Population 
    et_au
        Electric field
    vt_au
        Velocity 
    '''
    print(folder, version)
    if version == "new":
        if os.path.exists(folder + "/Time.txt"):
            t_au = np.loadtxt(folder + "/Time.txt", unpack=True)
        else:
            raise FileNotFoundError("Time file not found.")
        if os.path.exists(folder + "/Laser.txt"):
            Et_au = np.loadtxt(folder + "/Laser.txt", unpack=True)
        if os.path.exists(folder + "/Population.txt"):
            Pt = np.loadtxt(folder + "/Population.txt", unpack=True)
        if os.path.exists(folder + "/Velocity.txt"):
            Vtxr, Vtxi, Vtyr, Vtyi, Vtzr, Vtzi = np.loadtxt(folder + "/Velocity.txt", unpack=True)
            #Vt_au = 1j*np.zeros((3, t_au.shape[0]))
            Vt_au = np.array([Vtxr+1j*Vtxi, Vtyr+1j*Vtyi, Vtzr+1j*Vtzi])
    
    if version == "old":
        if os.path.exists(folder + "/Losses.txt"):
            Pt = np.loadtxt(folder + "/Losses.txt", unpack=True)
            t_au = Pt[0,:]/(constants.physical_constants["atomic unit of time"][0]*1.e+15)
            Pt = Pt[1:, :]
        else:
            print("Losses.txt does not exist")
        if os.path.exists(folder + "/EF.txt"):
            Et_au = np.loadtxt(folder + "/EF.txt", unpack=True)
            Et_au = Et_au[1:, :]
        else:
            print("/EF.txt does not exist")
        if os.path.exists(folder + "/J1.txt"):
            _, J1xr, J1xi, J1yr, J1yi, J1zr, J1zi = np.loadtxt(folder + "/J1.txt", unpack=True)
            J1 = np.array([J1xr+1j*J1xi, J1yr+1j*J1yi, J1zr+1j*J1zi])
        else:
            print("/J1.txt does not exist")
        if os.path.exists(folder + "/J2.txt"):
            _, J2x, J2y, J2z = np.loadtxt(folder + "/J2.txt", unpack=True)
            J2 = np.array([J2x, J2y, J2z])
        else:
            print("/J2.txt does not exist")
            #Vt_au = 1j*np.zeros((3, t_au.shape[0]))
        Vt_au = (J1+J2)

    #be sure everything has same dimensions -mainly if you analyse data while running a simulation
    dim   = np.min((Pt.shape[1], t_au.shape[0], Et_au.shape[1], Vt_au.shape[1]))
    Et_au = Et_au[:,:dim]
    t_au  = t_au[:dim]
    Pt    = Pt[:,:dim]
    Vt_au = Vt_au[:,:dim]
    if print_info:
        print(t_au.shape, Pt.shape, Vt_au.shape, Et_au.shape)
    return t_au, Pt, Et_au, Vt_au


import h5py
import json

def read_blockmatrix(filename, node, label, comm_size, dims):
    f = h5py.File(filename,'r+')[str(node)][label]
    H = []
    for key in range(comm_size):
        print(key)
        H.append(np.array(f[str(key)]['local']))
    H= np.array(H)
    H.dtype = complex
    H = np.reshape( H, dims ) 
    return H

def read_recap(filename):
    f = h5py.File(filename,'r+')
    js = f["json"]
    print(js.keys())
    dict_ = {}
    dict_["tb_file"]           = js["tb_file"][0]
    dict_["coulomb"]           = js["coulomb"][0]
    dict_["grid"]              = np.array([(js["grid"][0]), (js["grid"][1]), (js["grid"][2])])
    dict_["filledbands"]       = js["filledbands"][0]
    print(dict_["filledbands"])
    dict_["dt"]                = js["dt"][0]
    dict_["printresolution"]   = js["printresolution"][0]
    dict_["initialtime"]       = js["initialtime"][0]
    dict_["finaltime"]         = js["finaltime"][0]
    dict_["num_bands"]         = js["num_bands"][0]
    dict_["num_kpoints"]       = js["num_kpoints"][0]
    dict_["solver"]            = js["solver"][0]
    dict_["comm_size"]         = js["comm_size"][0]
    dict_["kpoints"]           = np.reshape( np.array(js["kpoints"]), (dict_["num_kpoints"],3) )
    dict_["A"]                 = np.reshape( np.array(js["A"]), (3,3) )
    dict_["B"]                 = np.reshape( np.array(js["B"]), (3,3) )

    H0=[]
    DM=[]
#    H0 = read_blockmatrix( filename, "H0", dict_["comm_size"], (dict_["num_kpoints"], dict_["num_bands"], dict_["num_bands"]) )
#    DM = read_blockmatrix( filename, "DM", dict_["comm_size"], (dict_["num_kpoints"], dict_["num_bands"], dict_["num_bands"]) )
    return dict_, H0, DM




def read_json(filename):
    with open(filename, 'r') as file:
        data = json.load(file)
    return data
