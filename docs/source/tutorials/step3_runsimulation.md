# Running a simulation
After you manage to install the code and get an input file, you can run a simulation by simply typing in the terminal: 
```bash 
/path/to/EDUS_2.0/build/EDUS <name>.json
```
During execution, there is a preprocessing part in which the input is read and all the EDUS variables are initialized accordingly. All the initialized variables are then written in the output file. If something looks weird, you can make sure that all the variables are initialized correctly in the code. <br>
After initialization, a folder called `Output` will be created in the folder where you run the simulation and inside it, you can find:
- `Time.txt`: the time grid in a.u.
- `Population.txt`: The number of holes in valence band and of electrons in conduction bands, normalized with respect to the number of k points of the simulation
- `Velocity.txt`: The expected value of the velocity operator over the state of the simulation: 

  $$
  \langle {\bf v} \rangle (t) = Tr(v \rho(t)) 
  $$

  There are six columns, that are <br>
  $Re(v_x)$, $Im(v_x)$, $Re(v_y)$, $Im(v_y)$, $Re(v_z)$, $Im(v_z)$
  <br>
  Being an expected value, if the imaginary part gets important something is going on in the simulation. Also, notice that the velocity should be zero at equilibrium.
- `output.h5`: This file contains all the information of the simulation, and can be later analysed using the postprocessing scripts

## Tutorial

### hBN IPA and HSEX spectra


Using the txt files generated, even during execution you can analyse your data using the python script that you can find in `/path/to/EDUS_2.0/PostProces/Recap.py` running it: 
```bash 
python3 /path/to/EDUS_2.0/PostProces/Recap.py --folder /path/to/Output/
```
