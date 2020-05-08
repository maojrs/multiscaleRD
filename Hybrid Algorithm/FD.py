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

'''
This code returns the solution of the reaction-diffusion equation laplace u=D*u_t+r*u with homogeneous 
Neumann boundary conditions and the initial condition at time 0 u_0:=u(x,0).
The solution is calculated witht he Finite Difference scheme and presents two different versions
of iteration matrices. To obtain the final concentration in each cell, we average over the 4 neighbouring
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
D=diffusion coefficient
r=first order macroscopic rate
L=domain length
'''

''' Possible functions as initial condition u_0'''

def u1(x):
    
    return exp(-0.5*((x[0]-12)**2+(x[1]-5)**2))*10

def u1sin(x): 
    
    return ((np.sin(x[0]+x[1])+1))
def u2(x):
	if x[0]<6:
		ux=0
	else:
		ux=exp(-(1/(0.1))*((x[0]-7)**2+(x[1]-6)**2))*(1/np.sqrt(0.1))
		
	return ux

def fconstant(x):
	if 6.5<=x[0]<=8.5 and 5<=x[1]<=7 :
		ux=30
	else:
		ux=0
	
	return ux    

def fconstant2(x):
	if 5<x[0]<6 and 2<x[1]<3 :
		ux=5
	else:
		ux=0
	
	return ux  

'Parameters'
u_0=fconstant

l=10
m=2*l-1  
m=l 
deltat=0.01 
n=1000 
L=12 
D=0.5 
r1=0.1 

h=L/(l-1) # grid-size


'Laplace Matrix with Neumann boundary conditions everywhere'

# Version 1

C=np.zeros((l, l))

for i in range(l):
    for j in range(l):
        if i==j:
            C[i,j]=4
        if abs(i-j)==1:
            
            C[i,j]=-1
           
C[0, 1]=-2
C[l-1, l-2]=-2            



A1 = scipy.sparse.bmat([[C if i == j  else -2*np.identity(l) if abs(i-j)==1
                        and j==0 and i==1 else -2*np.identity(l) if abs(i-j)==1 and j==l-1 else -np.identity(l) 
                        if abs(i-j)==1 else None for i in range(m)] for j in range(m)], format='bsr').toarray()
# Version 2
    
C2=np.zeros((l, l))

for i in range(l):
    for j in range(l):
        if i==j:
            C2[i,j]=4
        if abs(i-j)==1:
            
            C2[i,j]=-1
            
C2[0, 0]=2
C2[l-1, l-1]=2
J=np.zeros((l, l))

for i in range(l):
    for j in range(l):
        if i==j:
            J[i,j]=1
           
J[0, 0]=1/2
J[l-1, l-1]=1/2

A2 = scipy.sparse.bmat([[ C2 if i > 0 and i < l-1 and i==j else C2*0.5 if j==0 and i==0 else C2*0.5 if i==l-1 and j==l-1
                        else -J if abs(i-j)==1 else None for i in range(l)] for j in range(l)], format='bsr').toarray()

''' Create Solution Vector (column-wise)'''
   
x=np.linspace(0,L,m+1)
y=np.linspace(0,L,l+1)

U0=np.zeros(int((l)*(m))) # solution vector at time 0
a=0
for i in range(m):
    for j in range(l):
        U0[a]=u_0(np.array([x[i], y[j]])) # Initial Condition at time 0
        a=a+1

listU=[] # list of solution vectors for each time-step
listU.append(U0)

'''Strang-Splitting with implicit Euler'''

U=U0
A=A2 # choose a version
B=inv(np.identity(int(l*m))+D*deltat/(h**2)*A)  #iteration matrix
for t in range(n-1):
    U=U+deltat*0.5*r1*U
    U=B.dot(U)
    U=U+deltat*0.5*r1*U

    listU.append(U)
    
''' Translate solution vector into solution matrix'''

listM=[] # list of matrices
for t in range(n): 
    if t%1==0: # which time-steps we want to save
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
Prey=listM
np.save('FDSolution.npy', Prey) 

