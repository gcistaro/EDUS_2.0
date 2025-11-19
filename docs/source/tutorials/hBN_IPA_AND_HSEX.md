# hBN HSEX
You can find the input file in `/path/to/EDUS_2.0/test/hBN_HSEX.json`. <br>
In the input, notice that: 
- we excite the system with a laser intensity that give access to linear regime ($10^8 W/cm^2$), the laser is short (there is only 1 cycle and the period is 0.69 fs, being the frequency 6 eV, meaning that the zero-to-zero energy span is around 12eV).
- we use `dt=0.01fs` so that in the cycles of the laser we resolve with enough points, around 70, that is enough in linear regime
- we use a small grid of $20^2$ points as we are in linear regime
- the laser start when the simulation starts, being `t0 = 0 fs` and `initialtime = 0 fs`. 
- the tb_model is a 2D Tight-Binding model with an hexagonal lattice, with a gap of 7.25eV and a hopping parameter between neighbor atoms of -2.3eV, and only one orbital per atomic site (there are two atomic sites). This describes quite well the pz orbitals of hBN. 
- for the Coulomb interaction, we use a Rytova-Keldysh potential with $\epsilon=1.0$ and $r_0=10 ang$. 
<br>
Send the simulation running: 

```bash 
/path/to/EDUS_2.0/build/EDUS hBN_HSEX.json
```

## Analysing recap file

Using the txt files generated, even during execution you can analyse your data using the python script that you can find in `/path/to/EDUS_2.0/PostProces/Recap.py` running it: 

```bash 
python3 /path/to/EDUS_2.0/PostProces/Recap.py --folder /path/to/Output/
```

Here is how it should look like.

```{image} images/recap.png
:width: 100%
:align: center
```

Left and right columns are the same plot in time space and frequency space. <br>
In the rows, from top to bottom, there is the population in time for each band, the electric field in time and the mean value of the velocity over the state. <br>
What we can learn: 
- The electric field acts in the beginning of the simulations and creates a non-equilibrium population, electrons jump from valence to conduction leaving holes. After the pulse finish, the coherent population that is created is not decaying 
- From the velocity operator plot in time, you can see that there are oscillations with the frequency of the electric field, but also differenet frequencies contribute in the oscillations. This is reflected in the Fourier transform on the right, when you can clearly see some atomic-like excitations (bound excitons).
<br>

## Get the absorbance
The absorbance can be directly calculated from the electric field and the velocity operator using the formula:

$$
    \tilde{\sigma}(\omega) = 2 \frac{\omega}{\epsilon_0 c} \text{Im} \Biggl\{ \frac{\tilde{\boldsymbol{\varepsilon}}(\omega) \cdot \tilde{\textbf{d}}(\omega)}{\left| \tilde{\boldsymbol{\varepsilon}}(\omega) \right|^2} \Biggl\},
$$

This is what the script `Absorbance.py`. You can launch it using:

```bash
python3 /path/to/EDUS_2.0/PostProces/Absorbance.py --smearing=0.2 --folder=/path/to/Output/
```

The smearing is essential to make the velocity operator decay in time, so that the fft is well defined and the ripples present in the Fourier Transform obtained during the previous step will disappear. The output will look like this: 

```{image} images/absorbance.png
:width: 100%
:align: center
```

We can notice that:
- There is a strong bound exciton at energy 5.35 eV with a lot of oscillator strength (it is the dominant peak).
- There are other bound excitons that are not well resolved (to resolve them, we should send the simulation further in time and then decrease the smearing)



