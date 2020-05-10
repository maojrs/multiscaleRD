# -*- coding: utf-8 -*-
"""
Created on Mon May  4 11:42:00 2020

@author: Margarita
"""

from __future__ import division
import numpy as np
from math import exp

def movement(Particles,  deltat, D, Lx, Ly):
    
    ''' 
    Returns the new position of particles as a list after one time-step of length deltat:
    Particles=list of 2D arrays
    deltat=time-step size
    D= diffusion coefficient
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


def proliferation(Particles, rate, deltat):
    
    ''' 
    Simulates the first oder reaction A->2A or the proliferation of a species. It returns a list of positions (2D array)
    of new particles (children).
    Particles=list of 2D arrays
    rate=microscopic reaction rate of the first order reaction
    deltat=time-step size
    '''
  
    pproliferation=1-np.exp(-rate*deltat) 

    children=[] 
    for r in range(len(Particles)): 
        if pproliferation > np.random.rand():
           children.append(np.array(Particles[r]))
    return children

def euclid(v1, v2):
    
    '''
    Returns the euclidean distance between arrays v1 and v2 in an
    efficient way
    v1,v2=vectors with the same length
    '''
    dist = [(a-b)**2 for a,b in zip(v1, v2)]
    dist = np.sqrt(sum(dist))
    return dist

def dying(Particles,rate, deltat, NotImmune):
    
    '''
    Simulates the dying of one species, i.e. the reaction  A-> zero. However only
    species in the list 'NotImmune' can die.
    It returns a new list of arrays (position of the species) with removed
    dead species.
    Particles=list of 2D arrays (positions)
    rate=microscopic reaction rate
    deltat=time-step size
    NotImmune= positition of Particles that actually are able to die
    If all particles can die, just set NotImmune=Particles
    '''
    
    pdying=1-np.exp(-rate*deltat) 
    for r in reversed(range(len(Particles))):
        if pdying > np.random.rand() and any((Particles[r] == x).all() for x in NotImmune):
            Particles.pop(r)
    return Particles


def second_order_reaction(A, B, rate, deltat, sigma):
    
    '''Simulates the second order reaction A+B->2B, if A and B are closer then sigma 
    (reaction radius). It returns the new list of A particles with removed particles, a list of 
    B particles with new particles and a list of the new particles (new list of B=previous list of B +children)
    A, B=list of 2D arrays
    rate=microscopic reaction rate
    deltat=time-step size
    sigma=reaction radi
    '''
    p=1-exp(-deltat*rate)
    children=[]
    for i in reversed(range(len(A))):
        for j in reversed(range(len(B))):
            if euclid(A[i], B[j])<= sigma and p>np.random.rand():
                children.append(A[i])
                B.pop(j)
                A.pop(i)
                break
                 
    return A, children, B

def virtual(L, deltar, N,i):
    
    '''
    Assigns to virtual particles at the boundary in cells a position,
    such that they can react with each other in the function 'eatcompact'. Returns a list of 2D arrays.
    L=x-boundary coordinate
    deltar=length of the boundary
    N=number of particles we assign a position to
    i=boundary cell
    '''
    
    Virtual=[None]*N
    for j in range(N):
        Virtual[j]=np.array([np.random.uniform(L-deltar, L), np.random.uniform(deltar*i,deltar*(i+1))])
    return Virtual

def eatcompact(A, B, L, deltar, BC1, BC2, rate, sigma, deltat):
    
    '''
    Hybrid algorithm for second order reaction A+B->2B in and cross boundaries.
    It returns the new list of A particles with removed particles, a list of 
    B particles with new particles and a list of the new particles (new list of B=previous list of B +children)
    A, B=list of 2D arrays
    deltar=reaction radius
    BC1=boundary concentration of A particles
    BC2=boundary concentration of B particles
    rate=microscopic reaction rate
    sigma=reaction radius
    L=x-boundary coordinate
    deltat=time-step size
    '''

    for i in range(len(BC2)): # test for every  boundary cell (big cell)
        dec=(BC2[i]-int(BC2[i])) # crete virtual predators that can react with preys in the particle domain
        Virtual=virtual(L, deltar, int(dec),i)
        A, childrenbd1, B1=second_order_reaction(A, Virtual, rate, deltat, sigma)
        Virtual=virtual(L, deltar, 1,i)
        A, childrenbd2, B2=second_order_reaction(A, Virtual, rate*dec, deltat, sigma)
    A, children, B3=second_order_reaction(A, B, rate, deltat, sigma)
    
    kindergarten=childrenbd1+childrenbd2+children
    B=B1+B2+B3
    
    return A, kindergarten, B

