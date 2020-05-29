#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 15:46:15 2020

@author: bzfkostr
"""

import numpy as np

def fconstant(x):
	if 6.5<=x[0]<=8.5 and 5<=x[1]<=7 :
		ux=30
	else:
		ux=0
	
	return ux    


'Parameters'

'''
Space discretization
l+1=number of grid-cells in y direction
m+1=number of grid-cells in x direction

Time discretization
deltat=time-step size
timesteps=number of iterations 

Mathematical parameters
D=diffusion coefficient
r=first order macroscopic rate
a=domain length
l_coupling=number of cells in the FD scheme after the averaging
a=vetical length of the doman
L=horizonal length of the domain equals to the half of a
deltar=width of the boundary domain (one boundary cell has the legnth and width of deltar)

Computatinal Parameters
sim_number=number of simulations in each parallel simulation
numSimulations=how many simulations should run in parallel
s=which timestep we want to save
'''

u_0=fconstant
l=100+1 
m=l 
deltat=0.01 
timesteps=1000 
D=0.5 
r1=0.1 
a=12
  

''' Parameters only for Coupling'''
L=a/2 # half of the grid
deltar=np.sqrt(deltat*D*2) # boundary size
l_coupling=l-1 
sim_number=20 
numSimulations = 30
s=0.1 