import numpy as np


def cart(kcrys, B):
    #multiply A*kcrys with kcrys being (3 x num_k) array
    return np.matmul(B, kcrys)

def crys(kcart, invB):
    return np.matmul(invB, kcart)

def reduce(kcrys, limits):
    if( np.abs( limits[1]-limits[0]-1) > 1.e-06 ):
        print("Error in limits, they don't span a region of 1 :", limits[1]-limits[0])
        exit(1)
    
    #make sure everything is positive
    kcrys_reduced = kcrys + 1000

    #mod function to restrict in [0,1)
    kcrys_reduced = np.mod(kcrys, 1) + limits[0]
    return kcrys_reduced
    
