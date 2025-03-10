import numpy as np
from scipy.special import struve, yn

##import matplotlib.pyplot as plt
from scipy import constants
import sys
import os
import math


# read wannier
def read_wannierTB(filename):
    f = open(filename).readlines()
    fsplit = [line.split() for line in f]

    unit_cell = np.array(
        [
            [float(fsplit[1][0]), float(fsplit[1][1]), float(fsplit[1][2])],
            [float(fsplit[2][0]), float(fsplit[2][1]), float(fsplit[2][2])],
            [float(fsplit[3][0]), float(fsplit[3][1]), float(fsplit[3][2])],
        ]
    )
    ###print("\n Unit cell: ", unit_cell)
    num_bands = int(fsplit[4][0])
    ###print("\n num_bands: ", num_bands)
    num_R = int(fsplit[5][0])
    ###print("\n num_R: ", num_R)

    # skip lines with degeneracies
    index = 8 + int(num_R / 15)

    # read Hamiltonian
    R = np.zeros([num_R, 3])
    H = 1j * np.zeros([num_R, num_bands, num_bands])
    for iR in range(num_R):
        for ix in range(3):
            R[iR][ix] = fsplit[index][ix]
        index = index + 1
        for irow in range(num_bands):
            for icol in range(num_bands):
                H[iR][irow][icol] = float(fsplit[index][2]) + 1j * float(
                    fsplit[index][3]
                )
                index = index + 1
        index = index + 1

    # read Position
    r = 1j * np.zeros([3, num_R, num_bands, num_bands])
    for iR in range(num_R):
        index = index + 1
        for irow in range(num_bands):
            for icol in range(num_bands):
                r[0][iR][irow][icol] = float(fsplit[index][2]) + 1j * float(
                    fsplit[index][3]
                )
                r[1][iR][irow][icol] = float(fsplit[index][4]) + 1j * float(
                    fsplit[index][5]
                )
                r[2][iR][irow][icol] = float(fsplit[index][6]) + 1j * float(
                    fsplit[index][7]
                )
                index = index + 1

        index = index + 1
    return unit_cell, num_R, num_bands, R, H, r


# r0 = 10.
r0 = float(sys.argv[5])
###print(r0)
r0_au = r0 * 1.0e-10 / constants.physical_constants["Bohr radius"][0]
epsilon = float(sys.argv[6])
###print(epsilon)

# calculate the Rytova-Keldysh potential on the grid given
grid = [int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])]
name_tb = sys.argv[4]
A, num_R, num_bands, R, H, r = read_wannierTB(name_tb)
R_cart = np.matmul(R, A)
index_origin = np.where(np.linalg.norm(R_cart - np.zeros([1, 3]), axis=1) < 1.0e-07)
r_atom = np.array([r[0][index_origin], r[1][index_origin], r[2][index_origin]])
###print(r_atom.shape)
r_atom = np.reshape(r_atom, [3, num_bands, num_bands])
r_atom = np.array([np.diag(r_atom[0]), np.diag(r_atom[1]), np.diag(r_atom[2])])
r_atom = np.transpose(r_atom)


# gamma centered grid
R = np.array(np.meshgrid(range(grid[0]), range(grid[1]), range(grid[2])))
R = np.reshape(R, [3, grid[0] * grid[1] * grid[2]])
R = np.array([R[1, :], R[0, :], R[2, :]])
R = np.transpose(R)
###print(R.shape)

for iR, R_ in enumerate(R):
    for ix in range(3):
        if R[iR][ix] >= int(grid[ix] / 2) and grid[ix] > 1:
            R[iR][ix] = R[iR][ix] - int(grid[ix])
R = np.matmul(R, A)


###print("Calculating rytova keldysh")
RytovaKeldysh = np.zeros([R.shape[0], num_bands, num_bands])
for iR, R_ in enumerate(R):
    for ibnd1 in range(num_bands):
        for ibnd2 in range(num_bands):
            rnorm = np.linalg.norm(r_atom[ibnd1] - R_ - r_atom[ibnd2])
            RytovaKeldysh[iR][ibnd1][ibnd2] = (
                np.pi / (epsilon * r0_au) * (struve(0, rnorm / r0) - yn(0, rnorm / r0))
            )
            if rnorm < 1.0e-05:
                RytovaKeldysh[iR][ibnd1][ibnd2] = 0.0
            if math.isnan(RytovaKeldysh[iR][ibnd1][ibnd2]):
                RytovaKeldysh[iR][ibnd1][ibnd2] = (
                    np.pi
                    / (epsilon * r0_au)
                    * (struve(0, (rnorm + 0.01) / r0) - yn(0, (rnorm + 0.01) / r0))
                )
                if math.isnan(RytovaKeldysh[iR][ibnd1][ibnd2]):
                    print(
                        "R= ",
                        R_,
                        "atom positions: ",
                        r_atom[ibnd1],
                        r_atom[ibnd2],
                        "rnorm : ",
                        rnorm,
                        struve(0, rnorm / r0),
                        yn(0, rnorm / r0),
                        "is nan",
                    )
###print("Done")
###print(np.max(RytovaKeldysh), np.min(RytovaKeldysh))

###print("Saving RK potential in " + os.getcwd() + "/RytovaKeldysh.txt")
f = open(os.getcwd() + "/RytovaKeldysh.txt", "w")
f.write(
    "#File generated from "
    + name_tb
    + "with a grid "
    + str(grid[0])
    + "x"
    + str(grid[1])
    + "x"
    + str(grid[2])
    + "\n"
)
for iR, R_ in enumerate(R):
    for ibnd1 in range(num_bands):
        for ibnd2 in range(num_bands):
            #== f.write("{:15.8f}".format(R_[0]))
            #== f.write("{:15.8f}".format(R_[1]))
            #== f.write("{:15.8f}".format(R_[2]))
            #== f.write("{:15.8f}".format(np.real(r_atom[ibnd1][0])))
            #== f.write("{:15.8f}".format(np.real(r_atom[ibnd1][1])))
            #== f.write("{:15.8f}".format(np.real(r_atom[ibnd1][2])))
            #== f.write("{:15.8f}".format(np.real(r_atom[ibnd2][0])))
            #== f.write("{:15.8f}".format(np.real(r_atom[ibnd2][1])))
            #== f.write("{:15.8f}".format(np.real(r_atom[ibnd2][2])))
            f.write("{:15.8f}".format(RytovaKeldysh[iR][ibnd1][ibnd2]))
            f.write("\n")
f.close()
