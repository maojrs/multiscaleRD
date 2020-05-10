# -*- coding: utf-8 -*-
"""
Created on Mon May  4 12:12:19 2020

@author: Margarita

The following functions are made to inject particles from a conctinuous domain into the particle-based domain with rate gamma.
It requires the boundaryconcentration of the boundary of width deltar and length L (the domain length).
The functions will be used in the Strang splitting in the main file (Reaction+Injection, Diffusion, Reaction+Injection)
To guarantee an accurate result, we let the particles in the boundary alo proliferate (virual proliferation).
However this particles are not allowed to be injected in the same (Strang splitting) time-step (concentrationmovement0). On the other hand the particles
that had proliferated from the PREVIOUS time-step can be injected (concentrationmovement1).

"""
import numpy as np



def concentrationmovement0( Boundaryconcentration_t, deltat,deltar, L, Extra, gamma): 
    
    '''
    Returns a list of positions (2D arrays) of the new particles in the PBS domain (Children). The particles are injected with proabibility gamma. 
    Only the particles that did not proliferated in the same time-step can be injected. Therefore we have to subtract the 'Extra' particles 
    from the total boundaryconcentrations.
    deltat=time
    Booundaryconcentration_t=list or array of Boundaryconcentration of each cell of length and width deltar (L/deltar= number of boundary cells).
    deltat=time-step size
    gamme=injection rate
    deltar=boundary cell length and width
    Extra=number of proliferated particles in the boundary cell
    L=domain size
    '''
    Children=[]
    Pr=1-np.exp(-gamma*deltat) # probability of injection
    for i in range(len(Boundaryconcentration_t)):       
        Boundaryconcentration_t[i]=Boundaryconcentration_t[i]-Extra[i] # don't inject particles that just have proliferated in the boundary cell        
        integ, dec = int(np.floor(Boundaryconcentration_t[i])), Boundaryconcentration_t[i]-int(np.floor(Boundaryconcentration_t[i]))
        for v in range(integ): # test for every molecule           
            if Pr > np.random.rand():
                Children.append(np.array([np.random.uniform(L-deltar, L), np.random.uniform(deltar*i,deltar*(i+1))]))        
        if 1-np.exp(-gamma*deltat*dec)>np.random.rand(): # for the 'half' particle
            Children.append(np.array([np.random.uniform(L-deltar,L), np.random.uniform(deltar*i, deltar*(i+1))]))

    return Children


def concentrationmovement1( Boundaryconcentration_t, deltat,deltar,L, Extra, New, gamma):
    
    '''
    Returns a list of positions (2D arrays) of the new particles in the PBS domain (Children). The particles are injected with proabibility gamma. 
    The particles that did not proliferated in the same time-step and the proliferated particles from the last time step
    can be injected. Therefore we have to subtract the 'Extra' particles (proliferated in this time-step) and add the 'New' particles (proliferated
    in the previous time-step.) from and to the the total boundaryconcentrations.
    deltat=time
    Booundaryconcentration_t=list or array of Boundaryconcentration of each cell of length and width deltar (L/deltar= number of boundary cells).
    deltat=time-step size
    gamme=injection rate
    deltar=boundary cell length and width
    Extra=number of proliferated particles in the boundary cell
    New=number of prolieferated particles from the LAST time-step in the boundary cell
    L=domain size
    '''
    Children=[]
    Pr=1-np.exp(-gamma*deltat)
    for i in range(len(Boundaryconcentration_t)):
        Boundaryconcentration_t[i]=Boundaryconcentration_t[i]+Extra[i]-New[i] # add new preys from previous virtual proliferation (Extra), but not from the last one (New) 
        integ, dec = int(np.floor(Boundaryconcentration_t[i])), Boundaryconcentration_t[i]-int(np.floor(Boundaryconcentration_t[i]))
        for v in range(integ): # test for every molecule
            if Pr > np.random.rand():
                Children.append(np.array([np.random.uniform(L-deltar, L), np.random.uniform(deltar*i,deltar*(i+1))]))
        if 1-np.exp(-gamma*deltat*dec)>np.random.rand(): # for the 'half' particle
            Children.append(np.array([np.random.uniform(L-deltar,L), np.random.uniform(deltar*i, deltar*(i+1))]))

    return Children

'''Proliferation of virtual particles'''

def virtualproliferation(Boundaryconcentration_t, deltat, proliferationrate): 
    '''
    Returns the number of proliferated particles in the boundary cell for each boundary cell as an array.
    Booundaryconcentration_t=list or array of Boundaryconcentration of each cell of length and width deltar (L/deltar= number of boundary cells).
    deltat=time-step size
    proliferationrate=microscopic reaction rate
    '''
    virtualchildren=np.zeros(len(Boundaryconcentration_t))
    pproliferation=1-np.exp(-proliferationrate*deltat)  # probability of proliferation
    for i in range(len(Boundaryconcentration_t)): # test for every  boundary cell 
        integ, dec = int(np.floor(Boundaryconcentration_t[i])), Boundaryconcentration_t[i]-int(np.floor(Boundaryconcentration_t[i]))
        for r in range(integ): 
            if pproliferation > np.random.rand():
                virtualchildren[i]=virtualchildren[i]+1
        if 1-np.exp(-proliferationrate*deltat*dec) >np.random.rand():
            virtualchildren[i]=virtualchildren[i]+1
            
    return virtualchildren
