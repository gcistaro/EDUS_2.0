# Installation 

## Downloading the code 
EDUS is open source and can be downloaded from the github repository:
```bash
# Download EDUS
git clone https://github.com/gcistaro/EDUS_2.0.git
```

## Cmake installation
EDUS can be installed using the CMakeLists.txt file in the main directory. <br>
From a terminal, type:
```bash
# Install EDUS
mkdir build
cd build
cmake ..
make
```
Cmake will look for the dependencies of EDUS. In particular, make sure you have 
- A gcc compiler (for example, g++ or icpc)
- fftw installed 
- Intel-mkl installed
- mpi installed
- hdf5 library

If it sounds complicated to get all this before EDUS installation, you can refer to the next subsection to automatically install everything.

## Spack installation
To install EDUS with SPACK, first of all you need a version of spack in your computer

```bash
git clone -c feature.manyFiles=true https://github.com/spack/spack.git ~/spack
cd spack
git checkout v0.23.1
```

And need to source the setup file in spack:

```bash
. ~/spack/share/spack/setup-env.sh
```

Now spack is available. Try to find some packages that are usually already installed

```bash
spack external find autoconf automake git libtool m4 perl tar
spack compiler add
```

Now that you have spack, you can create a new environment. For example, create a folder:

```bash
cd
mkdir envs
cd envs
mkdir EDUS
cd EDUS
```

This will be the folder with your environment. Inside it, create a file called spack.yaml and inside copy this:

```bash
# This is a Spack Environment file.
#
# It describes a set of packages to be installed, along with
# configuration settings.
spack:
  # add package specs to the `specs` list
  specs:
  - EDUS@1.0%gcc@12.3.0 build_type=Debug
  view: true
  concretizer:
    unify: true
  develop:
    EDUS:
      path: path/to/EDUS_2.0
      spec: EDUS@1.0
  repos:
  - path/to/EDUS_2.0/spack/
```

Be careful to change path/to/ with the folder where the code can be found. Activate the environment you have just created and concretize it:

```bash
spacktivate . -p
spack concretize 
```

After this is done, you can install.

```bash
spack install -v 
```

