# Input variables
Information about the input variables can be found in `/path/to/EDUS_2.0/src/InputVariables/input_schema.json` where each variable is defined and described. Here we classify the variables depending on what they do.

## Light - matter variables
These variables serve to define the system that we are simulating.
- For the electronic system, we can identify: 
  - `"tb_model"`: the seedname of the tb_file, you can also give a path to the file
  - `"coulomb"`: the level of theory to use in the simulation. If true, we do HSEX; if false, we do IPA.
  - `"epsilon"` and `"r0"`: The parameters for the Rytova-Keldysh potential, in case `coulomb` is true
  - `"opengap"`: allows to open a gap (or eventually close it), with the number that you give
  - `"filledbands"`: number of bands to be filled with electrons. Important!! Make sure you fill fully a manifold of bands, or the wannierization idea will not work and you will soon get divergences in the time dynamics

- For the laser, we identify the following variables in an array defined inside a node called `"laser"`:
  - `"intensity"`: the peak intensity of the laser
  - `"frequency"` or `"wavelength"`: The frequency or wavelength of the laser, in the units defined
  - `"cycles"`: The number of cycles of the envelope, defined as a $sin^2$
  - `"polarization"`: the polarization of the laser in cartesian coordinates. If the norm is not 1, it will be normalized during preprocessing
  - `"t0"`: starting time of the laser, useful if you want to do pump-probe experiments
  - `"phase"`: phase of the laser. Put 0. if you want a sin pulse, $\frac{\pi}{2}$ if you want a cos pulse
## Convergence variables
These variables can be used to converge the numerical solver of the differential equation.
- `"grid"`: the k/R grid to be used in the time-dynamics, this is in general finer than the one used in DFT and everything is obtained using Wannier interpolation
- `"dt"`: the time step to use in the time evolution. In general, it is a good idea to put many time points (order of 100) in a laser cycle.
<br>
Depending on the intensity of the laser, you may need finer time resolution and more k points in the grid. <br>
For linear response (meaning that the intensity is lower than $10^9 W/cm^2$) we expect a grid of $16$ k points in each direction is enough to be converged, and `dt` can be such that you resolve a period with around 40 points in the time grid:

$$
dt = \frac{T}{40} = \frac{2\pi}{40\omega}
$$

## Simulation variables
- `"initialtime"` and `"finaltime"`: when the simulation starts and end in time

## Plot variables
- `"printresolution"`: Number of steps after we compute the observables and print them in txt and hdf5 files.
- `"printresolution_pulse`": Mainly to avoid to waste too much memory in the hdf5. When this number is lower than `printresolution`, the code will use this value for printing the big arrays in the hdf5 only when the laser is on, so we expect more dynamics to happen.
- `"toprint"`: array in which you can define: 
  - `"DMk_wannier"`: print in the hdf5 the density matrix in the wannier gauge
  - `"DMk_bloch"`: print in the hdf5 the density matrix in the eigenstate gauge
  - `"fullH"`: print in the hdf5 $H = H_0 + E(t)\cdot \xi + H_{ee}(t)$
  - `"SelfEnergy"`: print in the hdf5 $H_{ee}(t)$  
The big arrays can be then read and analysed using the python scripts or the mathematica script.
- `"kpath"`: an array of k points in crystal coordinates. If the array is defined, the band structure will be plotted in the path defined by the points in the array.