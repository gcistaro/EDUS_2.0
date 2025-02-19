from custom_functions.read import read_blockmatrix, read_recap
from custom_functions.basis import reduce
from custom_functions.basis import crys

import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt

dic, _, _ = read_recap("output.h5")


#create square in k space where we print
kpoints = reduce(dic["kpoints"],(0,1))
k=np.arange(-4,4,kpoints[1,1]-kpoints[0,1])
#k=np.arange(-1,1,kpoints[1,1]-kpoints[0,1])
y,x = np.meshgrid(k, k) # creates 2 np array 
kcart_ = np.transpose(np.array([y.flatten(), x.flatten(), np.zeros(y.flatten().shape[0])])) # reshuffle meshgrid correctly
kcrys = np.transpose(crys(np.transpose(kcart_), inv(dic["B"])))
#kcrys=kcart_
kcrys_reduced = reduce(kcrys, (0,1))
print("kcrys_reduced", kcrys_reduced)
print("kpoints", kpoints)
index = np.array([np.where( (np.abs(kc - kpoints) < (k[1]-k[0]) ).all(axis=1) )[0][0] for kc in kcrys_reduced ])


filename = "output.h5"
node = -1
label = "Xk"#"Epsilonk"
comm_size = dic["comm_size"]
dims = (dic["num_kpoints"], dic["num_bands"], dic["num_bands"])
X0k = read_blockmatrix(filename, node, label, comm_size, dims)

X0k = X0k[index,:,:]
X0k= X0k.reshape((len(k),len(k),dic["num_bands"], dic["num_bands"]))

#create colormap
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.colors as colors

N = 256; rgb=255.;
color=[(1.,1.,1.), (1.,0.,0),(1.,127./rgb,0.),(1.,1.,0.),(0.,1.,0.),(0.,0.,1.),(148./rgb,0.,211./rgb),(0.,0.,0.) ]
N2=round(N/(len(color)-1))
vals = np.ones((N, 4))
for index in np.linspace(1,len(color)-1,len(color)-1):
    index=int(index); #print(index, color[index-1], color[index])
    vals[(index-1)*N2:np.min([(index)*N2,vals.shape[0]]), 0] = np.linspace(color[index-1][0], color[index][0], np.min([N2,vals.shape[0]-(index-1)*N2]))
    vals[(index-1)*N2:np.min([(index)*N2,vals.shape[0]]), 1] = np.linspace(color[index-1][1], color[index][1], np.min([N2,vals.shape[0]-(index-1)*N2]))
    vals[(index-1)*N2:np.min([(index)*N2,vals.shape[0]]), 2] = np.linspace(color[index-1][2], color[index][2], np.min([N2,vals.shape[0]-(index-1)*N2]))
wrb = ListedColormap(vals)
#############################

for iband in range(2):
    for jband in range(2):
        fig,ax= plt.subplots()
        cax = fig.add_axes([0.9,0.0,0.01,1.0])
        im = ax.scatter(y, x,  c=np.abs(X0k[:,:,iband,jband]), cmap=wrb)
        fig.colorbar(im, cax, 'vertical')
        ax.set_aspect('equal')
        plt.show()

for iband in range(2):
    for jband in range(2):
        fig,ax= plt.subplots()
        index0 = np.where(np.abs(k) < 1.e-05)
        #ax.scatter(k, X0k[index0, np.linspace(0,k.shape[0],num=k.shape[0], endpoint=False).astype(np.int32), iband, jband])
        ax.scatter(np.linalg.norm([y,x],axis=0), X0k[:, :, iband, jband])
        plt.show()
