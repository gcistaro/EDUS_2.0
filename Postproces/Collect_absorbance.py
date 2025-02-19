import numpy as np
import matplotlib.pyplot as plt
f = ["/scratch/gcistaro/test_EDUS_versions/hBN/hBN_coulomb_old_100/r010_pol100_1cycle_10eV_int1e6_dt0.01_100x100grid/Output",
    "/scratch/gcistaro/Calculations/hBN_2bands/80"]#["20", "40", "60", "80"]


for folder in f:
    freq_eV, Absorbance = np.loadtxt(folder+"/absorbance.txt", unpack = True)
    plt.plot(freq_eV, Absorbance, label=folder)
plt.legend()
plt.show()
