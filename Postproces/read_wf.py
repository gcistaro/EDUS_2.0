from custom_functions.read import read_blockmatrix, read_recap
from custom_functions.basis import reduce
from custom_functions.basis import crys

import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt

dic, _, _ = read_recap("output.h5")


#create square in k space where we print
kpoints = reduce(dic["kpoints"],(0,1))
k=np.arange(-1,1,kpoints[1,1]-kpoints[0,1])
#k=np.arange(-1,1,kpoints[1,1]-kpoints[0,1])
y,x = np.meshgrid(k, k) # creates 2 np array 
kcart_ = np.transpose(np.array([y.flatten(), x.flatten(), np.zeros(y.flatten().shape[0])])) # reshuffle meshgrid correctly
kcrys = np.transpose(crys(np.transpose(kcart_), inv(dic["B"])))
#kcrys=kcart_
kcrys_reduced = reduce(kcrys, (0,1))
print("kcrys_reduced", kcrys_reduced)
print("kpoints", kpoints)
index = np.array([np.where( (np.abs(kc - kpoints) < (k[1]-k[0]) ).all(axis=1) )[0][0] for kc in kcrys_reduced ])
print(index.shape)

filename = "output.h5"
node = -1
label = "Avc"#"Epsilonk"
comm_size = dic["comm_size"]
dims = (dic["num_kpoints"], dic["num_kpoints"])#dic["num_bands"], dic["num_bands"])
X0k = read_blockmatrix(filename, node, label, comm_size, dims)

X0k = X0k[index,:]
X0k= X0k.reshape((len(k),len(k),dic["num_kpoints"]))#,len(k),len(k),len(k)))#dic["num_bands"], dic["num_bands"]))

#create colormap
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.colors as colors

N = 256; rgb=255.;
color=[(1.,1.,1.), (1.,0.,0),(1.,127./rgb,0.),(1.,1.,0.),(0.,1.,0.),(0.,0.,1.),(148./rgb,0.,211./rgb),(0.,0.,0.) ]
N2=round(N/(len(color)-1))
vals = np.ones((N, 4))
for index_ in np.linspace(1,len(color)-1,len(color)-1):
    index_=int(index_); #print(index, color[index-1], color[index])
    vals[(index_-1)*N2:np.min([(index_)*N2,vals.shape[0]]), 0] = np.linspace(color[index_-1][0], color[index_][0], np.min([N2,vals.shape[0]-(index_-1)*N2]))
    vals[(index_-1)*N2:np.min([(index_)*N2,vals.shape[0]]), 1] = np.linspace(color[index_-1][1], color[index_][1], np.min([N2,vals.shape[0]-(index_-1)*N2]))
    vals[(index_-1)*N2:np.min([(index_)*N2,vals.shape[0]]), 2] = np.linspace(color[index_-1][2], color[index_][2], np.min([N2,vals.shape[0]-(index_-1)*N2]))
wrb = ListedColormap(vals)
#############################

fig,ax= plt.subplots()
iwf = 0
cax = fig.add_axes([0.9,0.0,0.01,1.0])
im = ax.scatter(y,x,  c=np.abs(X0k[:,:,iwf]), cmap=wrb)
fig.colorbar(im, cax, 'vertical')
ax.set_aspect('equal')
plt.show()


phase = np.zeros(index.shape[0])
print(phase.shape)
X0 = np.reshape(X0k[:,:,iwf], (index.shape[0]))

for ix0, x0_ in enumerate(X0):
    if np.abs(X0[ix0]) > 1.e-06:
        phase[ix0] = np.angle(x0_)
    else:
        phase[ix0] = 0.
phase = np.reshape(phase, (len(k),len(k)))

fig,ax= plt.subplots()
cax = fig.add_axes([0.9,0.0,0.01,1.0])
im = ax.scatter(y,x,  c=phase, cmap=wrb)
fig.colorbar(im, cax, 'vertical')
ax.set_aspect('equal')
plt.show()




##fig,ax= plt.subplots()
##cax = fig.add_axes([0.9,0.0,0.01,1.0])
##im = ax.scatter(y,x,  c=np.real(X0k[:,:,iwf]), cmap="gist_heat")
##fig.colorbar(im, cax, 'vertical')
##ax.set_aspect('equal')
##plt.show()
##fig,ax= plt.subplots()
##cax = fig.add_axes([0.9,0.0,0.01,1.0])
##im = ax.scatter(y,x,  c=np.imag(X0k[:,:,iwf]), cmap="gist_heat")
##fig.colorbar(im, cax, 'vertical')
##ax.set_aspect('equal')
##plt.show()