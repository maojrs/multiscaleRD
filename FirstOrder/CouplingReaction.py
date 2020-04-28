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


'''Parameters'''

D=0.5 #Diffusion coefficient
timesteps=1000-1   #Number of timesteps without zero
deltat=0.01       # time step size
l=100 #number of grid cells
a=12  #length of domain
L=a/2  #len of hybdrid domain
r1=0.1 # proliferation constant
x0 = np.array([L+1, L])  #location of source
deltar=np.sqrt(deltat*D*2) # length of cell edge at boundary layer for hybrid scheme (squares of dx by dx)
dx_hist =a/l # length of histogram cell edge (squares of dx_hist by dx_hist)
dtshould = deltar*deltar/(2.0*D) # timestep size
print(dtshould, deltat, 'Should be equal')

'''Calculate number of particles at the boundary cel'''

maxtime = deltat*(timesteps) # maximum time simulation can reach
Time=np.linspace(0, maxtime, timesteps+1)
listC=np.load('./1501F.npy') # gets Data from continuous solution
yarray = np.arange(0,a,deltar) # Array to locate boundary cells  
averageNumberParticles = np.zeros((len(yarray),timesteps+1))

xlimits = [a/2,a/2+deltar]
for i in range(len(yarray)):

    ylimits = [yarray[i], yarray[i] + deltar]
    for k in range(timesteps+1):
        if k == 0:
            averageNumberParticles[i,k] =  0.0

        else:

            averageNumberParticles[i,k] = (deltar) * (deltar) * listC[k][int(i/(len(yarray)/l)),int(l/2) ] # X=C*V


Boundaryconcentration=averageNumberParticles

'''Diffusion'''

def movement( List,  deltat, D):
    Q=len(List)
    for i in reversed(range(Q)):
        List[i]=np.array(List[i])+np.random.normal(0,np.sqrt(2*deltat*D),2) # Euler-Maruyama scheme
        '''Boundary Conditions'''
        if List[i][0]>= L: 
            List.pop(i)
        else: # Reflecting
            if List[i][0]<=0:
                List[i][0]=-List[i][0]
                
            if List[i][1]>=a:
                List[i][1]=a-(List[i][1]-a)
            if List[i][1]<=0:
                List[i][1]=-List[i][1]
        
    return List

'''First Injection'''

def concentrationmovementfunction0( Boundaryconcentrationt, deltat,deltar, Extra):  # Extra= proliferated particles from bundary cell
    Children=[] # list of new particles from the injection
    gamma=D/((deltar)**2) # injection rate
    Pr=(1-np.exp(-gamma*deltat)) 
    for i in range(len(Boundaryconcentrationt)): 
        
        Boundaryconcentrationt[i]=Boundaryconcentrationt[i]-Extra[i] # don't inject particles that just have proliferated in the boundary cell
        
        integ, dec = int(np.floor(Boundaryconcentrationt[i])), Boundaryconcentrationt[i]-int(np.floor(Boundaryconcentrationt[i]))
        for v in range(integ): # test for every molecule
           
            if Pr > np.random.rand():
                Children.append(np.array([np.random.uniform(L-deltar, L), np.random.uniform(deltar*i,deltar*(i+1))]))
        
        if 1-np.exp(-gamma*deltat*dec)>np.random.rand(): # for the 'half' particle
            Children.append(np.array([np.random.uniform(L-deltar,L), np.random.uniform(deltar*i, deltar*(i+1))]))

    return Children

'''Second Injection'''

def concentrationmovementfunction1( Boundaryconcentrationt, deltat,deltar, Extra, New):
    Children=[]
    gamma=D/((deltar)**2)
    Pr=(1-np.exp(-gamma*deltat)) 
    for i in range(len(Boundaryconcentrationt)):
        
        Boundaryconcentrationt[i]=Boundaryconcentrationt[i]+Extra[i]-New[i] # add new preys from previous virtual proliferation (Extra), but not from the last one (New) 

        integ, dec = int(np.floor(Boundaryconcentrationt[i])), Boundaryconcentrationt[i]-int(np.floor(Boundaryconcentrationt[i]))
        
        for v in range(integ): # test for every molecule
           
            if Pr > np.random.rand():
                Children.append(np.array([np.random.uniform(L-deltar, L), np.random.uniform(deltar*i,deltar*(i+1))]))
        
        if 1-np.exp(-gamma*deltat*dec)>np.random.rand(): # for the 'half' particle
            Children.append(np.array([np.random.uniform(L-deltar,L), np.random.uniform(deltar*i, deltar*(i+1))]))

    return Children

'''Proliferation of virtual particles'''

def virtualproliferation(Boundaryconcentrationt, deltat,deltar, proliferationrate): 
    virtualchildrenlist=np.zeros(len(Boundaryconcentrationt))
    pproliferation=1-np.exp(-proliferationrate*deltat) 
    for i in range(len(Boundaryconcentrationt)): # test for every  boundary cell 
        integ, dec = int(np.floor(Boundaryconcentrationt[i])), Boundaryconcentrationt[i]-int(np.floor(Boundaryconcentrationt[i]))
        virtualchildren=[]
        for r in range(integ): 
            if pproliferation > np.random.rand():
                virtualchildren.append(np.array([np.random.uniform(L,L+deltar), np.random.uniform(deltar*i, deltar*(i+1))])) # asign position is actually unnecessary
                
        
        if 1-np.exp(-proliferationrate*deltat*dec) >np.random.rand():
            virtualchildren.append(np.array([np.random.uniform(L,L+deltar), np.random.uniform(deltar*i, deltar*(i+1))])) # asign position is actually unnecessary
        virtualchildrenlist[i]=len(virtualchildren)

    return virtualchildrenlist

'''First Order Reaction'''

def proliferation(listPrey, proliferationrate, deltat):
    pproliferation=1-np.exp(-proliferationrate*deltat) 

    listchildren=[] # Position list of new molecules from thr reaction
    for r in range(len(listPrey)): 
        if pproliferation > np.random.rand():
            listchildren.append(np.array(listPrey[r]))
    return listchildren

''' Doubled Stang-Splitting'''

def functionsimulation(simulations):

    PreySimulation=simulations*[None]
    for s in range(simulations):
    	PreyPosition=[] # list of particles (2d arrays)
    	PreyPositionHalfTime=[] # list of particles at each saved time step
    	Reference=[] #list of reference solutions at each time step
    	for t in range(timesteps+1):

    		'''Injection'''
    		PreyChildrenVirtual=virtualproliferation(Boundaryconcentration[:,t], deltat,deltar, r1)
    		PreyChildrenI1= concentrationmovementfunction0(Boundaryconcentration[:,t], deltat*1/2,deltar, PreyChildrenVirtual)

    		'''Reaction'''
    		PreyChildrenR1= proliferation(PreyPosition, r1, deltat*1/2) # Reaction step

    		'''Diffusion'''
    		PreyPosition = PreyPosition + PreyChildrenI1 + PreyChildrenR1 # add injected and proliferated preys to the population
    		PreyPosition = movement(PreyPosition,  deltat, D)
    		'''Reaction'''
    		PreyChildrenR2 = proliferation(PreyPosition, r1, deltat*1/2)

    		'''Injection'''
    		New=virtualproliferation(Boundaryconcentration[:,t], deltat,deltar, r1)
    		PreyChildrenI2 = concentrationmovementfunction1(Boundaryconcentration[:,t], deltat*1/2,deltar, PreyChildrenVirtual, New)

    		PreyPosition= PreyPosition + PreyChildrenR2 + PreyChildrenI2 # add injected and proliferated preys to the population

    		if t*deltat in np.arange(0, maxtime+deltat, 0.1 ): 
        		PreyPositionHalfTime.append(PreyPosition)

        		if s==simulations-1:
        			Reference.append(listC[t]) 

    	print( s, len(PreyPosition),t)

    	PreySimulation[s]=PreyPositionHalfTime
        
    return PreySimulation, Reference

'''Multi-Processing'''   

def runParallelSims(simnumber):

    # Define seed
    np.random.seed(simnumber)
   

    # Save to file
    PreySimulation, Reference=functionsimulation(100)
    np.save( '/nfs/datanumerik/bzfkostr/Codes/dataNew/ReactionParticles'+str(simnumber)+'.npy', PreySimulation)
    np.save( '/nfs/datanumerik/bzfkostr/Codes/dataNew/ReactionReference'+str(simnumber)+'.npy', Reference)
	
	# run simulation
    print("Simulation " + str(simnumber) + ", done.")

numSimulations = 30 #30
num_cores = multiprocessing.cpu_count()
print(num_cores), 'Kerne'
pool = Pool(processes=num_cores)
iterator = [i for i in range(numSimulations)]
pool.map(partial(runParallelSims), iterator)

