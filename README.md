# MultiscaleRD
Code to accompany the paper ''Coupling particle-based reaction-diffusion 
simulations with reservoirs mediated by reaction-diffusion PDEs.'' 
by M. Kostré, C. Schütte, F. Noé and M. J. del Razo.
In this code we implemented a hybrid algorithm that couples particle-based simulations with spatially and time dependent reservoir modeled with reaction-diffusion PDEs.
It calculates the solution of the 2D reaction-diffusion equation (RDE) with a mix of the Finite Difference and the Crank-Nicolson algorithm, runs particle based reaction-diffusion simulations, and couples both codes to each other by injecting particles
from the concentration reservoir - given by the solution of the PDE, into the particle-based domain, where it runs a particle-based simulation. The examples of this code show the consistency between the micro and macroscopic level 
up to the second order. 
We implemented two different scripts:
* one for reactions up to the first order (Coupling.py)
* one for reactions up to the second order (CouplingLV.py)
The last named simulates the classical Lotka-Volterra dynamics. 

## Prerequisites

* python
* git: to clone this repository 

## Example

This example shows the reference and hybrid solution for the preys in the Lotka-Volterra dynamics.

![Watch the video](Videos/PreyReferenceVideo.gif)
![Watch the video](Videos/PreyHybridVideo.gif)


## How to run this code?

* Clone this repository to your local machine: git clone https://github.com/MargKos/multiscaleRD.git
* Solve the RDE by running FD.py, this gives the reference solution(s)
* Run many simulation with Coupling.py or CouplingLV.py
* Calculate the mean-field of the PBS and create plots with Plot.py


## Folder Organization

* In the FD.py file we implemented the Finite Difference scheme for a RDE with homogeneous Neumann boundary conditions. 
* The Coupling.py file runs many (here 30) parallel simulations. We recommend to adjust this file in such a way that the number
of parallel simulations corresponds to the number of cores of the computer. 
* Injection.py and Reaction.py content the functions used in Coupling.py: Functions from Injection.py return lists of positions of new particles injected from the 
continuous domain. These lists are appended to the list of particles in Coupling.py. In Reaction.py function for a  particle-based simulations (PBS) up to the second order are implemented, like proliferation (A->2A), dying (A->0) or the 
classical movement of particles calculated from the Euler-Muruyama scheme. The functions are used in Coupling.py and returns for example list of positions of newly created particles.
* Plot.py calculates the mean-field creates hybrid plots, where the right part corresponds to the reference solution obtained through
FD.py and the left part to the average of the trajectories densities obtained from Coupling.py. The discretization parameters needs to be specified there according to the setting of FD.



## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
