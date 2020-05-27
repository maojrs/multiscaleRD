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


'''
This main file couples a PBS to a mean-field concentration and returns trajectories of particles for different simulations. 
Before running this file it is necessary to calculate the continuous solution with FD.py. The functions for reactions are from Reaction.py and 
the one for injection from Injection.py
The input is:
D1, D2=Diffusion coefficient of A and B
timesteps=stepsnumber of time
deltat=time-step size
l=number of cells in the FD scheme
a=vetical length of the doman
L=horizonal length of the domain
r1=first order microscopic rate A-> 2A
r2_macro=second order macroscopic rate, that corresponds to the rate in the FD scheme
r2=second order microscopic rate A+B-> 2B
sigma=reaction radius
r3=zero reaction B->0 miroscopic rate
deltar=width of the boundary domain (one boundary cell has the legnth and width of deltar)
The code consits of the following components
1) Calculate the boundary concentration for EACH timestep
2) Iteration with Strang Splitting using the function from Reaction.py and Injection.py: Injection, Reaction, Diffusion, Reaction, Injection.
3) Multiprocessing that does many simulations at same time
'''


'''Parameters'''

timesteps=4000
deltat=0.0025      
l=100
a=10
L=a/2  
D1=0.3
D2=0.1
r1=0.15 
sigma=0.01
r2_macro=0.05
r2=r2_macro/(np.pi*(sigma**2)) 
r3=0.1

'''1. Calculate the boundary councentration 'Boundaryconcentration' from the FD solution '''


x0 = np.array([L+1, L])  #location of source
deltar1=np.sqrt(deltat*D1*2) # length of cell edge at boundary layer for hybrid scheme (squares of dx by dx) for species A
deltar2=np.sqrt(deltat*D2*2) # length of cell edge at boundary layer for hybrid scheme (squares of dx by dx) for species B
dx_hist =a/l # length of histogram cell edge (squares of dx_hist by dx_hist)
yarray1 = np.arange(0,a,deltar1) # Array to lcoate boundary cells for species A
yarray2 = np.arange(0,a,deltar2) # Array to lcoate boundary cells for species B 

# Simulation parameters
dtshould1 = deltar1*deltar1/(2.0*D1) # time-step size
dtshould2 = deltar2*deltar2/(2.0*D2) # time-step size
print(dtshould1, dtshould2, deltat, 'Should be equal')
maxtime = deltat*(timesteps) # maximum time simulation can reach
Time=np.linspace(0, maxtime, timesteps+1)

listC1=np.load('./Data/FDSolution1.npy') # gets Data from continuous solution of A
listC2=np.load('./Data/FDSolution2.npy') # gets Data from continuous solution of B

averageNumberParticles1 = np.zeros((len(yarray1),timesteps+1))
averageNumberParticles2 = np.zeros((len(yarray2),timesteps+1))
xlimits1 = [a/2,a/2+deltar1]
xlimits2 = [a/2,a/2+deltar2]

gamma1=D1/((deltar1)**2)  #injection rate
gamma2=D2/((deltar2)**2)

'''Boundary concentration for A'''

for i in range(len(yarray1)):

    ylimits1 = [yarray1[i], yarray1[i] + deltar1]

    for k in range(timesteps+1):
        if k == 0:
            averageNumberParticles1[i,k] =  0.0

        else:

            averageNumberParticles1[i,k] = (deltar1) * (deltar1) * listC1[k][int(i/(len(yarray1)/l)),int(l/2) ]


'''Boundary concentration for B'''

for i in range(len(yarray2)):

    ylimits2 = [yarray2[i], yarray2[i] + deltar2]
    for k in range(timesteps+1):
        if k == 0:

            averageNumberParticles2[i,k] =  0.0
        else:

            averageNumberParticles2[i,k] = (deltar2) * (deltar2) * listC2[k][int(i/(len(yarray2)/l)),int(l/2) ]


Boundaryconcentration1=averageNumberParticles1
Boundaryconcentration2=averageNumberParticles2

''' 2. Iteration'''

def functionsimulation(simulations, ts):
    '''
    Returns each lists for A and B consisting of trajectories at every desirable time-step and for every simulation and the reference solutions
    simulations=number of simulations
    ts=which timestep we save: 0, ts, ts*2
    
    PreySimulation=saves all simulations
    PreyPostion=updating list containing the CURRENT position of preys as a 
    2D array
    PreyPositionHalfTime=contains the positions at each desirable time step
    Structure: PreySimulation=[[Simulation_1], [Simulation_2]...], Simulation_1=[[Time_1], [Time_2],..], for example Time_1=[[1,2], [3,4],[1.4,4]] (positions)
    Analogou for predator
    '''

    PreySimulation=simulations*[None]
    PredatorSimulation=simulations*[None]
    for s in range(simulations):
    	PreyPosition=[]
    	PreyPositionHalfTime=[]
    	PredatorPosition=[]
    	PredatorPositionHalfTime=[]
    	Reference1=[]
    	Reference2=[]
    	for t in range(timesteps+1): #Prey has a position
    		

    		'''Injection'''
    		PreyChildrenInj= concentrationmovement0(Boundaryconcentration1[:,t], deltat*1/2,deltar1,L,np.zeros(len(Boundaryconcentration1[:,t])), gamma1)
    		PredChildrenInj= concentrationmovement0(Boundaryconcentration2[:,t], deltat*1/2,deltar2,L,np.zeros(len(Boundaryconcentration2[:,t])), gamma2)

    		'''Reaction'''
    		PreyChildrenProlif=proliferation(PreyPosition, r1, deltat*0.25)
    		PredatorPosition=dying(PredatorPosition,r3, deltat*0.25, PredatorPosition)
    		PreyPosition, PredChildrenReact, PredB=eatcompact(PreyPosition, PredatorPosition, L, deltat*0.5,Boundaryconcentration1[:,t],Boundaryconcentration2[:,t], r2, sigma, deltat)
    		PreyChildrenProlif=proliferation(PreyPosition, r1, deltat*0.25)
    		PredatorPosition=dying(PredatorPosition,r3, deltat*0.25, PredB)

    		'''Put them all together'''
    		PreyPosition=PreyChildrenInj+PreyChildrenProlif+PreyPosition
    		PredatorPosition=PredChildrenInj+PredatorPosition+PredChildrenReact
            
    		'''Diffusion'''
    		PreyPosition=movement(PreyPosition,  deltat, D1, L, a)
    		PredatorPosition=movement(PredatorPosition,  deltat, D2, L, a)
        
    		"Reaction"
    		PreyChildrenProlif=proliferation(PreyPosition, r1, deltat*0.25)
    		PredatorPosition=dying(PredatorPosition,r3, deltat*0.25, PredatorPosition)
    		PreyPosition, PredChildrenReact, PredB=eatcompact(PreyPosition, PredatorPosition, L, deltat*0.5,Boundaryconcentration1[:,t],Boundaryconcentration2[:,t], r2, sigma, deltat)
    		PreyChildrenProlif=proliferation(PreyPosition, r1, deltat*0.25)
    		PredatorPosition=dying(PredatorPosition,r3, deltat*0.25, PredB)

    		PreyChildrenInj= concentrationmovement0(Boundaryconcentration1[:,t], deltat*1/2,deltar1,L,np.zeros(len(Boundaryconcentration1[:,t])), gamma1)
    		PredChildrenInj= concentrationmovement0(Boundaryconcentration2[:,t], deltat*1/2,deltar2,L,np.zeros(len(Boundaryconcentration2[:,t])), gamma2)
    		'''Put them all together'''
    		PreyPosition=PreyChildrenInj+PreyChildrenProlif+PreyPosition
    		PredatorPosition=PredChildrenInj+PredatorPosition+PredChildrenReact
    		if np.round(t*deltat,15) in np.round(np.arange(0, maxtime, ts),15):
        		PreyPositionHalfTime.append(PreyPosition)
        		PredatorPositionHalfTime.append(PredatorPosition)
        		if s==simulations-1:
        			Reference1.append(listC1[t])
        			Reference2.append(listC2[t])
        		print(t*deltat)


    	
    	print( s, len(PreyPosition), len(PredatorPosition),t)

    	PreySimulation[s]=PreyPositionHalfTime
    	PredatorSimulation[s]=PredatorPositionHalfTime
    return PreySimulation, Reference1, PredatorSimulation, Reference2


'''3. Multi-Processing'''   

sim_number=100
PreySimulation, Reference1, PredatorSimulation, Reference2=functionsimulation(1, 0.1)
np.save( './Simulation/Reference1.npy', Reference1) # saves reference solution at the CORRECT time-step
np.save( './Simulation/Reference2.npy', Reference2) # saves reference solution at the CORRECT time-step

def runParallelSims(simnumber):

    # Define seed
    np.random.seed(simnumber)

    # Save to file
    PreySimulation, Reference1, PredatorSimulation, Reference2=functionsimulation(sim_number, 0.1)
    np.save( './Simulation/PreyParticles'+str(simnumber)+'.npy', PreySimulation)
    np.save( './Simulation/PredatorParticles'+str(simnumber)+'.npy', PredatorSimulation)
  
	# run simulation
    print("Simulation " + str(simnumber) + ", done.")

numSimulations = 30 # how many simulations should run in parallel
num_cores = multiprocessing.cpu_count()
print(num_cores), 'Kerne'
pool = Pool(processes=num_cores)
iterator = [i for i in range(numSimulations)]
pool.map(partial(runParallelSims), iterator)

