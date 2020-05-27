# MultiscaleRD
Code to accompany the paper ''Coupling particle-based reaction-diffusion 
simulations with reservoirs mediated by reaction-diffusion PDEs.'' 
by M. Kostré, C. Schütte, F. Noé and M. J. del Razo.

In this code we implemented a hybrid scheme that couples particle-based reaction-diffusion simulations with spatially and time dependent reservoirs mediated by reaction-diffusion PDEs. It solves the reaction-diffusion PDE with a finite difference sheme with a Crank-Nicolson integrator, and it implements particle based reaction-diffusion simulations based on the Doi model, similar as how it is done in [ReaDDy2](https://readdy.github.io/). The hybrid scheme consistently couples the particle-based simulationthe to the PDE-mediated reservoir. We verify the scheme using two different examples: 
* A diffusion process with a proliferation reaction. This can be generalized for any systems with zeroth- and/or first-order reactions. (`Hybrid Algorithm/FD.py` and `Hybrid Algorithm/Coupling.py`)
* A Lotka-Volterra (predator-prey) reaction diffusion process. This can be generalized to any system with up to second-order reactions. (`Hybrid Algorithm/FD_LV.py` and `Hybrid Algorithm/CouplingLV.py`)

## Requirements

* python 3.x
* Numpy 1.12
* Scipy 0.19
* matplotlib 2.0 or above
* Multiprocessing 0.7
* git: to clone this repository 

## How to run this code?

* Clone this repository to your local machine: git clone https://github.com/MargKos/multiscaleRD.git
* Solve the reaction-diffusion PDE of each example by running the finite difference code, this gives the reference solution(s)
* Run many particle-based simulations with Coupling.py or CouplingLV.py
* Calculate the mean-field of the PBS and create plots with Plot.py

## Folder Organization

* In `Hybrid Algorithm/FD.py` we implemented the Finite Difference scheme for a RDE with homogeneous Neumann boundary conditions. 
* The `Hybrid Algorithm/Coupling.py` file runs many (here 30) parallel simulations. We recommend to adjust this file in such a way that the number
of parallel simulations corresponds to the number of cores of the computer. 
* `Hybrid Algorithm/Injection.py` and `Hybrid Algorithm/Reaction.py` content the functions used in Coupling.py: Functions from Injection.py return lists of positions of new particles injected from the 
continuous domain. These lists are appended to the list of particles in Coupling.py. In Reaction.py function for a  particle-based simulations (PBS) up to the second order are implemented, like proliferation (A->2A), dying (A->0) or the 
classical movement of particles calculated from the Euler-Muruyama scheme. The functions are used in Coupling.py and returns for example list of positions of newly created particles.
* `Hybrid Algorithm/Plot.py` calculates the mean-field creates hybrid plots, where the right part corresponds to the reference solution obtained through `Hybrid Algorithm/FD.py` and the left part to the average of the trajectories densities obtained from `Hybrid Algorithm/Coupling.py`. The discretization parameters needs to be specified there according to the settings in `Hybrid Algorithm/FD.py`.

## Sample solution

This example shows videos of the reference and hybrid solution for the preys in the Lotka-Volterra system (predator-prey).

<img src="Videos/PreyReferenceVideo.gif" width="400"> <img src="Videos/PreyHybridVideo.gif" width="400" />

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
