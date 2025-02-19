# tb_file
name of the file _tb.dat
 - **type**: name of the file _tb.dat
 - **default**: name of the file _tb.dat
 - **title**: name of the file _tb.dat

# fermienergy
Fermi energy of the system
 - **type**: Fermi energy of the system
 - **default**: Fermi energy of the system
 - **title**: Fermi energy of the system

# fermienergy_units
Units used for the Fermi energy of the system
 - **type**: Units used for the Fermi energy of the system
 - **default**: Units used for the Fermi energy of the system
 - **enum**: Units used for the Fermi energy of the system
 - **title**: Units used for the Fermi energy of the system

# coulomb
If true, we go beyond IPA using HSEX theory
 - **type**: If true, we go beyond IPA using HSEX theory
 - **default**: If true, we go beyond IPA using HSEX theory
 - **title**: If true, we go beyond IPA using HSEX theory

# grid
Dimension of the k grid (R as well) used for the simulation
 - **type**: Dimension of the k grid (R as well) used for the simulation
 - **items**: Dimension of the k grid (R as well) used for the simulation
 - **minItems**: Dimension of the k grid (R as well) used for the simulation
 - **maxItems**: Dimension of the k grid (R as well) used for the simulation
 - **default**: Dimension of the k grid (R as well) used for the simulation
 - **title**: Dimension of the k grid (R as well) used for the simulation

# filledbands
Number of filled bands for the simulation (ideally will be either this of fermienergy
 - **type**: Number of filled bands for the simulation (ideally will be either this of fermienergy
 - **default**: Number of filled bands for the simulation (ideally will be either this of fermienergy
 - **title**: Number of filled bands for the simulation (ideally will be either this of fermienergy

# dt
Time resolution used for propagating Density Matrix
 - **type**: Time resolution used for propagating Density Matrix
 - **default**: Time resolution used for propagating Density Matrix
 - **title**: Time resolution used for propagating Density Matrix

# dt_units
Units for the Time resolution used for propagating Density Matrix
 - **type**: Units for the Time resolution used for propagating Density Matrix
 - **default**: Units for the Time resolution used for propagating Density Matrix
 - **enum**: Units for the Time resolution used for propagating Density Matrix
 - **title**: Units for the Time resolution used for propagating Density Matrix

# printresolution
Time steps needed to print in hdf5 and observables
 - **type**: Time steps needed to print in hdf5 and observables
 - **default**: Time steps needed to print in hdf5 and observables
 - **title**: Time steps needed to print in hdf5 and observables

# initialtime
Initial time of the simulation
 - **type**: Initial time of the simulation
 - **default**: Initial time of the simulation
 - **title**: Initial time of the simulation

# initialtime_units
Units for the Time resolution used for Initial time of the simulation
 - **type**: Units for the Time resolution used for Initial time of the simulation
 - **default**: Units for the Time resolution used for Initial time of the simulation
 - **enum**: Units for the Time resolution used for Initial time of the simulation
 - **title**: Units for the Time resolution used for Initial time of the simulation

# finaltime
Final time of the simulation
 - **type**: Final time of the simulation
 - **default**: Final time of the simulation
 - **title**: Final time of the simulation

# finaltime_units
Units for the Time resolution used for Final time of the simulation
 - **type**: Units for the Time resolution used for Final time of the simulation
 - **default**: Units for the Time resolution used for Final time of the simulation
 - **enum**: Units for the Time resolution used for Final time of the simulation
 - **title**: Units for the Time resolution used for Final time of the simulation

# solver
Solver used for time propagation of the Density Matrix
 - **type**: Solver used for time propagation of the Density Matrix
 - **default**: Solver used for time propagation of the Density Matrix
 - **enum**: Solver used for time propagation of the Density Matrix
 - **title**: Solver used for time propagation of the Density Matrix

# order
Order of the Solver used for time propagation of the Density Matrix
 - **type**: Order of the Solver used for time propagation of the Density Matrix
 - **default**: Order of the Solver used for time propagation of the Density Matrix
 - **title**: Order of the Solver used for time propagation of the Density Matrix

# lasers
description of all the lasers
 - **type**: description of all the lasers
 - **title**: description of all the lasers
 - **minItems**: description of all the lasers
 - **maxItems**: description of all the lasers
 - **items**: description of all the lasers

