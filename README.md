[![Static Badge](https://img.shields.io/badge/DOI-10.1021%2Facs.jctc.2c00674-blue?style=flat&logo=DOI)
](https://doi.org/10.1021/acs.jctc.2c00674)

# EDUS - Electron Dynamics and Ultrafast Spectroscopy

This code can be used to propagate the electronic density matrix of an extended system up to HSEX level. It can simulate the interaction of the system with a (classical) electric field. \
The equations we solve are:<br />
   $`i{{\partial \rho(\textbf{k})}\over{\partial t}} = \Big[ H_0(k) +E(t)\cdot \xi(k) + H_{ee}(k), \rho(k)\Big] + i E(t)\cdot \nabla_k \rho(k)`$

The theory behind the code plus details about the implementation can be found in our paper [Theoretical Approach for Electron Dynamics and Ultrafast Spectroscopy (EDUS)](https://doi.org/10.1021/acs.jctc.2c00674). If you find our paper or the code useful, please consider citing us.

---

## Building instructions
#### Requirements
Basics:
```bash
gcc13
cmake
git
```
#### Link dependencies
    • MKL (for lapack/blas routines, eventually for fftw)
    • MPI (openmpi/mpich)
    • fftw (mpi version if you want to compile with mpi)
    • hdf5 (for ouputing results in an output.h5 data file)

It is important to note that the `mpi`, `fftw`, and `hdf5` must have matching version numbers. I.e.,
```bash
gcc/13.2.0
openmpi_gcc/5.0.2_gcc13.2.0
fftw/3.3.10_mpi5_gcc13
hdf5/1.14.3_gcc13_mpi5
```


#### Building
For all distros:
```bash
git clone https://github.com/gcistaro/EDUS_2.0.git
cd EDUS_2.0
mkdir build; cd build
cmake ..
make
ctest
```

---

## Input file
A sample input file can be found in `<repository>/ci-test/inputs`. For now you need to define all the variables, we will soon assign default values to many of the variables. One thing is that, currently, the code supports just reading the coulomb interaction from a file. You need to produce it using the script in `PostProces/RytovaKeldysh.py <nk1> <nk2> <nk3> tb_file`.

---

## Usage

Lorem ipsum dolor sit amet, consectetur adipiscing elit. Curabitur non metus nisl. Duis mattis sit. 

---

## Adding new code
1. Before adding something new, make sure you are not repeating something is already there.
2. To add new code, create a new class, unless you don’t just want to add some features (methods) in a class that already exists. Whenever you want to push something in the remote folder, do not do it in the main branch, but create a new branch with small changes, only related to what you want to add. When you are done with the changes, the first step is to make sure that you did not break old parts of the code, unless you did not find a bug. So, make everything and run ctest. If your tests are failing, consider what you changed in the “old” code that was working before. 
3. To create a new branch, give meaningful names, so that if it stays there for years we will remember what we were doing in that branch.
4. IMPORTANT:  push only the files you want to really to change. To see what you changed, you can use the terminal command git status. If something changed and you did not change it, you can run `git diff <file>`, that will compare your local file with the remote one. If you did some changes (even spaces) that are not needed, you can restore the remote file with `git restore <file>`.
5. After your changes are done, create a pull request and we will add them in main.

### HOW TO ADD A CLASS
If you want to add a class from scratch, you need to add to the code two different files:
1. A header file .hpp, that you will link to all the parts of the code that need your class using `#include <”file”.hpp>`. This class contains the declaration of your class together with the declaration of all the methods. “Declaration” means a line where you define return type, name, arguments but never their definition. This will be included in many parts of the code, so make sure it will be compiled only once with #ifdef.
2. A source file .cpp, which contains the actual code to be used. You need to add this file in the list of __SOURCES in CMakeLists.txt so it will be compiled as an object file with the code. In this file you need to write the definition of everything you just declared in your .hpp file. 
