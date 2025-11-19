# Introduction

EDUS (Electron Dynamics for Ultrafast Spectroscopy) is a C++ code that can propagate the electronic density matrix in condensed matter systems.<br>
It solves the time-dependent equations for the density matrix (Von Neumann Equation):

$$
i{{\partial \rho(\textbf{k})}\over{\partial t}} = \Big[ H_0(k) +E(t)\cdot \xi(k) + H_{ee}(k), \rho(k)\Big] + i E(t)\cdot \nabla_k \rho(k)
$$ 

where: 

- $\rho(k)$ is the electronic density matrix
- $E(t)$ is the electric field that perturbs the system
- $H_0(k)$ is the unperturbed Hamiltonian, that can be obtained from a DFT calculation or a Tight Binding model, and is saved in the format _tb.dat (see more details in wannier90 documentation to understand the file)
- $\xi(k)$ is the Berry connection, defined as the Fourier transform of the position operator [Blount1962], also saved in the file _tb.dat
- $H_{ee}(k)$ is the electron-electron coupling, defined using the Coulomb interaction. There are currently three approximation that can be used: 
  - IPA (Independent Particle Approximation): No electron coupling  
  
  $$
  H_{ee} = 0
  $$

  - RPA (Random Phase Approximation): Only Hartree potential 
  
  $$
  H_{ee} = V_{Hartree}
  $$

  - HSEX (Hartree Plus Screened Exchange): Hartree potential and screened Fock: 
  
  $$
  H_{ee} = V_{Hartree} + W_{Fock}
  $$ 
  
  <br>
It must be noted that $V_{Hartree}$ is obtained from the bare Coulomb interaction and $W_{Fock}$ from the screened Coulomb interaction. The two quantities can be obtained from an effective model in the code, or read ab-initio using the kcw interface (see section for more details).
The code works in the Wannier basis, defined as a rotation of the eigenstates that makes the Bloch wavefunctions smooth with respect to the variable k.   


:::{note}
You can insert notes, warnings, and tips using MyST.
:::

## What you will learn
- How to install the code
- How to write an input file
- How to run a simulation and look at your results