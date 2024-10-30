import sys
import numpy as np

file1 = sys.argv[1]
file2 = sys.argv[2]

arg1 = np.loadtxt(file1, unpack=True)
arg1 = arg1[1:].flatten()

arg2 = np.loadtxt(file2, unpack=True)
arg2 = arg2[1:].flatten()

print("comparing ", file1, " with ", file2)
passed = True
for i in range(arg1.shape[0]):
    if ( np.abs( arg1[i]-arg2[i] ) > 1.e-09  ):
        if ( np.logical_and ( np.abs(arg1[i]) < 1.e-09 , np.abs(arg2[i]) < 1.e-09 ) ):
            print(arg1[i], arg2[i], np.abs(arg1[i]-arg2[i]))
            sys.exit(1)
            passed = False

if(passed):
    print("... Passed!")
else:
    print("... Not Passed!")
