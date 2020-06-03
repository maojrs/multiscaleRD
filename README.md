# MultiscaleRD
Code to accompany the paper ''Coupling particle-based reaction-diffusion 
simulations with reservoirs mediated by reaction-diffusion PDEs.'' 
by M. Kostré, C. Schütte, F. Noé and M. J. del Razo.
[[arXiv]](https://arxiv.org/pdf/2006.00003.pdf)

In this code we implemented a hybrid scheme that couples particle-based reaction-diffusion simulations with spatially and time dependent reservoirs mediated by reaction-diffusion PDEs. It solves the reaction-diffusion PDE usimg a finite difference sheme with a Crank-Nicolson integrator, and it implements particle based reaction-diffusion simulations based on the Doi model, similar as how it is done in [ReaDDy2](https://readdy.github.io/). The hybrid scheme consistently couples the particle-based simulation to the to the PDE-mediated reservoir. We verify the scheme using two examples: 
* A diffusion process with a proliferation reaction. This can be generalized for any systems with zeroth- and/or first-order reactions. (`multiscaleRD/FD_proliferation.py` and `multiscaleRD/Coupling_proliferation.py`)
* A Lotka-Volterra (predator-prey) reaction diffusion process. This can be generalized to any system with up to second-order reactions. (`mutliscaleRD/FD_LV.py` and `multiscaleRD/Coupling_LV.py`)

## Requirements

* python 3.x
* Numpy 1.12
* Scipy 0.19
* matplotlib 2.0 or above
* Multiprocessing 0.7
* git: to clone this repository 

## How to run this code?

* Clone this repository to your local machine: git clone https://github.com/MargKos/multiscaleRD.git
* If desired, change the default mathematical and numerical parameters for the finte difference scheme and the particle based simulation in `multiscaleRD/Parameters_(example).py`
* Solve the reaction-diffusion PDE of each example by running the finite difference code (`multiscaleRD/FD_(example).py`), this yields the reference solution(s) and the reservoir.
* Run several particle-based simulations by running `multiscaleRD/Coupling_(example).py`
* Run `multiscaleRD/Plot_(example).py` to calculate the average over particle-based simulations and to create the plots. 

## Folder Organization

* `multiscaleRD/Parameters_(example).py` contains mathematical parameters (micro and macroscopic rates, diffusion coefficients), numerical parameters (timestepsize, gridsize, boundarysize...) and computational parameters (number of parallel simulations).
* In `multiscaleRD/FD_(example).py` we implemented the Finite Difference scheme for a RDE with homogeneous Neumann boundary conditions. 
* The `multiscaleRD/Coupling_(example).py` file runs many (here 30) parallel simulations. We recommend to adjust this file in such a way that the number of parallel simulations corresponds to the number of cores of the computer used. 
* `multiscaleRD/Injection.py` and `multiscaleRD/Reaction.py` contain the functions to implement the injection and reaction procedures used by the hybrid scheme in `multiscaleRD/Coupling_(example).py`:
  * Functions from `multiscaleRD/Injection.py` returns lists of positions of new particles injected from the 
continuous domain. These lists are appended to the list of particles in `multiscaleRD/Coupling_(example).py`.
  * In `multiscaleRD/Reaction.py`, functions for  particle-based reaction-diffusion simulations (up to the second order) are implemented, such as proliferation (A->2A), degradation (A->0) or the diffusion of particles modeled by the Euler-Maruyama scheme. The functions return for example list of positions of newly created particles.
* `multiscaleRD/Plot_(example).py` calculates the average ofver the particle-based simulations, and it creates the hybrid plots and animations, where the right part corresponds to the reference solution obtained 
through `multiscaleRD/FD_(example).py` and the left part to the average of the trajectories densities obtained from `multiscaleRD/Coupling_(example).py`, see example below.
* The folder 'multiscaleRD/Simulations' contains all the particle trajectories resulting from the hybrid scheme simulation from `multiscaleRD/Coupling_(example).py`.
* The folder'multiscaleRD/Solutions' contains solutions of the reaction-diffusion PDE and the average over particle-based simulations calculated from the simulations obtained by 'multiscaleRD/Simulations'. We use the files here for plotting or analysing the solution in `multiscaleRD/Plot_(example).py`.

### Notes on code notation
* In the LV example the solutions of the preys are denoted by 1, where the solutions of the predators are denoted by 2.
* In the prolieferation example we have only one species, so we don't have a specific ending for the simulations files, e.g. instead of 'PreyParticles', we just write 'Particles'.

## Sample solution

This example shows two videos for the reaction-diffusion dynamics of the concentration of preys in the Lotka-Volterra system (predator-prey). The first video shows the reference simulation obtained with a finite-difference scheme. The second video shows the results of the hybrid simulation. The left half corresponds to the average concentration over several particle-based simulations using our scheme, and the right half corresponds to the PDE-mediate reservoir.  

<img src="Videos/PreyReferenceVideo.gif" width="400"> <img src="Videos/PreyHybridVideo.gif" width="400" />

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
