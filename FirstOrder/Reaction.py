# -*- coding: utf-8 -*-
"""
Created on Mon May  4 11:42:00 2020

@author: Margarita
"""

from __future__ import division
import numpy as np

def movement(Particles,  deltat, D, Lx, Ly):
    
    ''' 
    Returns the new position of particles as a list after one time-step of length deltat:
    Particles=list of 2D arrays
    D= diffusion coefficient
    deltat=time-step size
    Lx=length of the horozontal boundary
    Ly=length of the vetical boundary
    We use the  Euler-Maruyama scheme, and REFLECTING boundary conditions.
    '''
    
    for i in reversed(range(len(Particles))):
        Particles[i]=np.array(Particles[i])+np.random.normal(0,np.sqrt(2*deltat*D),2) 
        if Particles[i][0]>= Lx: 
            Particles.pop(i)
        else: 
            if Particles[i][0]<=0:
               Particles[i][0]=-Particles[i][0]
                
            if Particles[i][1]>=Ly:
                Particles[i][1]=Ly-(Particles[i][1]-Ly)
            if Particles[i][1]<=0:
                Particles[i][1]=-Particles[i][1]
        
    return Particles


def proliferation(Particles, proliferationrate, deltat):
    
    ''' 
    Simulates the first oder reaction A->2A or the proliferation of a species. It returns a list of positions (2D array)
    of new particles (children).
    Particles=list of 2D arrays
    deltat=time-step size
    proliferationrate=microscopic reaction rate of the first order reaction'''
  
    pproliferation=1-np.exp(-proliferationrate*deltat) 

    children=[] 
    for r in range(len(Particles)): 
        if pproliferation > np.random.rand():
           children.append(np.array(Particles[r]))
    return children

