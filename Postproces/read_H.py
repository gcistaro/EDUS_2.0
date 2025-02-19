from custom_functions.read import read_blockmatrix, read_recap
from custom_functions.basis import reduce
from custom_functions.basis import crys

import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt

dic, H0, DM = read_recap("output.h5")
#filename = "SelfEnergy100.h5"

#SelfEnergy = read_blockmatrix(filename)
itfin = int( (dic["finaltime"]-dic["initialtime"])/float(dic["printresolution"]) )+2 
print("itfin=", itfin)

#create square in k space where we print
k=np.arange(-2,2,0.05)
y,x = np.meshgrid(k, k) # creates 2 np array 
kcart_ = np.transpose(np.array([y.flatten(), x.flatten(), np.zeros(y.flatten().shape[0])])) # reshuffle meshgrid correctly
kcrys = np.transpose(crys(np.transpose(kcart_), inv(dic["B"])))
kcrys_reduced = reduce(kcrys, (0,1))
kpoints = dic["kpoints"]
print("kcrys_reduced", kcrys_reduced)
print("kpoints", kpoints)
index = np.array([np.where( (np.abs(kc - kpoints) < 1.e-02 ).all(axis=1) )[0][0] for kc in kcrys_reduced ])

###DM = DM[index,:,:]
###DM= DM.reshape((len(k),len(k),dic["num_bands"], dic["num_bands"]))
####H0=H0[:-1,:-1,:,:]
###print(H0.shape)
###for iband in range(2):
###    for jband in range(2):
###        fig,ax= plt.subplots()
###        ax.scatter(y, x,  c=DM[:,:,iband,jband])
###        ax.set_aspect('equal')
###        plt.show()

for it in range(itfin):
    filename = str(it) + ".h5"
    label = "Ht"
    comm_size = dic["comm_size"]
    dims = (dic["num_kpoints"], dic["num_bands"], dic["num_bands"])
    Ht = read_blockmatrix(filename, label, comm_size, dims)

    print(Ht)
