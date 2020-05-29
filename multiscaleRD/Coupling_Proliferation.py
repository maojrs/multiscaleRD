#/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 10:23:43 2019

@author: bzfkostr
"""

from __future__ import division
import numpy as np
import matplotlib
matplotlib.use("agg")
import time
import multiprocessing
from multiprocessing import Pool
from functools import partial
from Reaction import *
from Injection import *
from Parameters_Prolieferation import *

'''
This main file couples a PBS to a mean-field concentration and returns trajectories of particles for different simulations. 
Before running this file it is necessary to calculate the continuous solution with FD.py. The functions for reactions are from Reaction.py and 
the one for injection from Injection.py
The input is:
D=Diffusion 
timesteps=stepsnumber of time
deltat=time-step size
l=number of cells in the FD scheme
a=vetical length of the doman
L=horizonal length of the domain
r1=first order microscopic rate
deltar=width of the boundary domain (one boundary cell has the legnth and width of deltar)
The code consits of the following components
1) Calculate the boundary concentration for EACH timestep
2) Iteration with Strang Splitting using the function from Reaction.py and Injection.py: Injection, Reaction, Diffusion, Reaction, Injection
3) Multiprocessing that does many simulations at same time
'''


deltar=np.sqrt(deltat*D*2) 

x0 = np.array([L+1, L])  #location of source
dx_hist =a/l_coupling # length of histogram cell edge (squares of dx_hist by dx_hist)
dtshould = deltar*deltar/(2.0*D) # timestep size
print(dtshould, deltat, 'Should be equal')
gamma=D/((deltar)**2) # injection rate

'''1. Calculate the boundary councentration 'Boundaryconcentration' from the FD solution '''

maxtime = deltat*(timesteps-1) # maximum time simulation can reach
Time=np.linspace(0, maxtime, timesteps)
listC=np.load('./Solutions/FDSolution.npy') # gets Data from continuous solution
yarray = np.arange(0,a,deltar) # Array to locate boundary cells  
averageNumberParticles = np.zeros((len(yarray),timesteps))

xlimits = [a/2,a/2+deltar]
for i in range(len(yarray)):
    ylimits = [yarray[i], yarray[i] + deltar]
    for k in range(timesteps):
        if k == 0:
            averageNumberParticles[i,k] =  0.0
        else:
            averageNumberParticles[i,k] = (deltar) * (deltar) * listC[k][int(i/(len(yarray)/l_coupling)),int(l_coupling/2) ] # X=C*V

Boundaryconcentration=averageNumberParticles

''' 2. Iteration'''

def functionsimulation(simulations, ts):
    '''
    Returns a list consisting of trajectories at every desirable time-step and for every simulation and the reference solution
    simulations=number of simulations
    ts=which timestep we save: 0, ts, ts*2
    
    PreySimulation=saves all simulations
    PreyPostion=updating list containing the CURRENT position of preys as a 
    2D array
    PreyPositionHalfTime=contains the positions at each desirable time step
    Structure: PreySimulation=[[Simulation_1], [Simulation_2]...], Simulation_1=[[Time_1], [Time_2],..], for example Time_1=[[1,2], [3,4],[1.4,4]] (positions)
    '''
    PreySimulation=simulations*[None]
    for s in range(simulations):
    	PreyPosition=[] 
    	PreyPositionHalfTime=[] # list of particles at each saved time step
    	Reference=[] #list of reference solutions at each time step
    	for t in range(timesteps):

    		'''Injection'''
    		PreyChildrenVirtual=virtualproliferation(Boundaryconcentration[:,t],  r1,deltat)
    		PreyChildrenI1= concentrationmovement0(Boundaryconcentration[:,t], deltat*1/2,deltar,L, PreyChildrenVirtual, gamma)

    		'''Reaction'''
    		PreyChildrenR1= proliferation(PreyPosition, r1, deltat*1/2) # Reaction step

    		'''Diffusion'''
    		PreyPosition = PreyPosition + PreyChildrenI1 + PreyChildrenR1 # add injected and proliferated preys to the population
    		PreyPosition = movement(PreyPosition,  deltat, D, L, a)
    		'''Reaction'''
    		PreyChildrenR2 = proliferation(PreyPosition, r1, deltat*1/2)

    		'''Injection'''
    		New=virtualproliferation(Boundaryconcentration[:,t],  r1,deltat)
    		PreyChildrenI2 = concentrationmovement1(Boundaryconcentration[:,t], deltat*1/2,deltar, L, PreyChildrenVirtual, New, gamma)

    		PreyPosition= PreyPosition + PreyChildrenR2 + PreyChildrenI2 # add injected and proliferated preys to the population

    		if np.round(t*deltat,10) in np.round(np.arange(0, maxtime, ts),10): 
        		PreyPositionHalfTime.append(PreyPosition)
        		if s==simulations-1:
        			Reference.append(listC[t]) 

    	print( s, len(PreyPosition),t)
    	PreySimulation[s]=PreyPositionHalfTime
        
    return PreySimulation, Reference

'''3. Multi-Processing'''   



PreySimulation, Reference=functionsimulation(sim_number, s)
np.save( './Solutions/Reference.npy', Reference) # saves reference solution at the CORRECT time-step
def runParallelSims(simnumber):

    # Define seed
    np.random.seed(simnumber)
   

    # Save to file
    PreySimulation, Reference=functionsimulation(sim_number, 0.1)
    np.save( './Simulations/Particles'+str(simnumber)+'.npy', PreySimulation)
  
	# run simulation
    print("Simulation " + str(simnumber) + ", done.")


num_cores = multiprocessing.cpu_count()
print(num_cores), 'Kerne'
pool = Pool(processes=num_cores)
iterator = [i for i in range(numSimulations)]
pool.map(partial(runParallelSims), iterator)

