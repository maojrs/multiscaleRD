# MultiscaleRD
Code to accompany the paper ''Coupling particle-based reaction-diffusion 
simulations with reservoirs mediated by reaction-diffusion PDEs.'' 
by M. Kostré, C. Schütte, F. Noé and M. J. del Razo.
In this code we implemented a hybrid algorithm that couples particle-based simulations with spatially and time dependent reservoirs modeled with reaction-diffusion PDEs.
It calculates the solution of the 2D reaction-diffusion equation (RDE) with the Finite Difference scheme, simulates reaction-diffusion particle-based  and couples them to each other by injecting particles
from the concentration reservoir, given by the solution of the PDE, into the particle-based domain. The examples of this code show the consistency between the micro and macroscopic level 
up to the second order. 
We implemented two different scripts:
* one for reactions up to the first order (Coupling.py)
* one for reactions up to the second order (CouplingLV.py)
The last named simulates the classical Lotka-Volterra dynamics. 

## Prerequisites

* python
* git: to clone this repository 

```
Give examples
```

## How to run this code?

* Clone this repository to your local machine:...
* Solve the RDE by running FD.py, this gives the reference solution(s)
* Run many simulation with Coupling.py or CouplingLV.py
* Create plots with Plot.py


## Folder Organization

* In the FD.py file we implemented the Finite Difference scheme for a RDE with homogeneous Neumann boundary conditions. 
* The Coupling.py file runs many (here 30) parallel simulations. We recommend to adjust this file in such a way that the number
of parallel simulations corresponds to the number of cores of the computer. 
* Injection.py and Reaction.py content the functions used in Coupling.py: Functions from Injection.py return lists of positions of new particles injected from the 
continuous domain. These lists are appended to the list of particles in Coupling.py. In Reaction.py function for a  particle-based simulations (PBS) up to the second order are implemented. 
* Plot.py creates hybrid plots, where the right part correpsonds to the reference solution obtained through
FD.py and the left part to the average of the trajectories densities obtained by Coupling.py.

### Example

![Hybrid plot for the first order reaction](Reaction9Hybrid.pdf)

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc
