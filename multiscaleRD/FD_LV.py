#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 10:40:31 2019

@author: bzfkostr
"""
from __future__ import division
import numpy as np
import scipy.sparse
from numpy.linalg import inv
from math import exp
from Parameters_LV import *

'''
This code returns the solution of the reaction-diffusion equation laplace u=D*u_t+r(u1,u2) with r( being the reaction function of the Lotka-Volerra equation)
, where u=(u1,u2) with homogeneous 
Neumann boundary conditions and the initial conditions at time 0 u1_0:=u1(x,0) and u2_0:=u2(x,0).
The solution is calculated witht he Finite Difference scheme and Strang Splitting. To obtain the final concentration in each cell, we average over the 4 neighbouring
grid-cells of FD. This requires only one more cell in the calculation. 
In the script we also give some possible initial condition u_0.
The skript requires:
    
Space discretization
l+1=number of grid-cells in y direction
m+1=number of grid-cells in x direction

Time discretization
deltat=time-step size
n=number of iterations 

Mathematical parameters
D1, D2=diffusion coefficients of u1 and u2
r1=first order macroscopic rate
r2_macro=second oder macroscopic rate
r3=zero reaction mmacroscopic rate
L=domain length

REMINDER: Lotka-Volterra equation
A->2A with rate r1
A+B->2B with rate r2_macro
B->0 with rate r3

A=Preys
B=Predators

'''



'Laplace Matrix with Neumann boundary conditions everywhere'

C=np.zeros((l, l))

for i in range(l):
    for j in range(l):
        if i==j:
            C[i,j]=4
        if abs(i-j)==1:
            
            C[i,j]=-1
           
C[0, 1]=-2
C[l-1, l-2]=-2            



A = scipy.sparse.bmat([[C if i == j  else -2*np.identity(l) if abs(i-j)==1
                        and j==0 and i==1 else -2*np.identity(l) if abs(i-j)==1 and j==l-1 else -np.identity(l) 
                        if abs(i-j)==1 else None for i in range(m)] for j in range(m)], format='bsr').toarray()

''' Create Solution Vector (column-wise) for preys and predators'''    
    
x=np.linspace(0,a,m+1)
y=np.linspace(0,a,l+1)


U0Prey=np.zeros(int((l)*(m)))
a=0
for i in range(m):
    for j in range(l):
        U0Prey[a]=uPrey(np.array([x[i], y[j]]))
        a=a+1

U0Pred=np.zeros(int((l)*(m)))
a=0
for i in range(m):
    for j in range(l):
        U0Pred[a]=uPred(np.array([x[i], y[j]]))
        a=a+1

listt=np.linspace(0, timesteps*deltat, timesteps)

U1=np.zeros(int(l*m))
listU1=[]
listU1.append(U0Prey)
U1=U0Prey

U2=np.zeros(int(l*m))
listU2=[]
listU2.append(U0Pred)
U2=U0Pred

'''Iteration matrices'''

B1=inv(np.identity(int(l*m))+D1*deltat/(h**2)*A)
B2=inv(np.identity(int(l*m))+D2*deltat/(h**2)*A)


'''Strang-Splitting with implicit Euler'''

for t in range(timesteps-1):
    U1=U1+deltat*0.5*(r1*U1-r2_macro*U1*U2) #half timestep
    U2=U2+deltat*0.5*(r2_macro*U1*U2-r3*U2)
    U1=B1.dot(U1) #full timestep
    U2=B2.dot(U2)
    U1=U1+deltat*0.5*(r1*U1-r2_macro*U1*U2)
    U2=U2+deltat*0.5*(r2_macro*U1*U2-r3*U2)

    listU1.append(U1)
    listU2.append(U2)
    
def functionMatrix(listU):
    listM=[]
    for t in range(timesteps): 
        if t%1==0: # which time-steps you want to save
            helpM=np.zeros((l, m))
            Ut=listU[t]
            k=0
           
            for j in range(m):
                for i in range(l):
                    helpM[i,j]=Ut[k]
                    k=k+1
            M=np.zeros((l-1,m-1))
            for i in range(l-1):
                for j in range(m-1):
                    M[i,j]=(helpM[i+1,j]+helpM[i,j]+helpM[i,j+1]+helpM[i+1,j+1])/4
             
            listM.append(M)
    return listM
Prey=functionMatrix(listU1)
Pred=functionMatrix(listU2)

np.save('./Solutions/FDSolution1', Prey)
np.save('./Solutions/FDSolution2', Pred)
