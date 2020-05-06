#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 12:11:01 2019

@author: bzfkostr
"""

from __future__ import division
import numpy as np
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns
#import animation

'''
Returns heatmap plots. Therefore, it first calculates the particles density and then 
averages over all simulations.
'''

'''Parameters'''

l=120 # nuber of cells
a=12 # length of the domain
D=0.5 # Diffusion constant
L=a/2

dx_hist=a/l
deltat=0.01 # Time-step size

'''Densities'''

def Discretization(a, discret, Particles): 
    
    '''
    Return the concentration of particles.
    a=domain length
    discret=discrtization parameter (number of cells)
    Particles=list of 2D arrays
    '''
    
    h=a/discret

    xPositions=[]
    yPositions=[]
    for i in range(len(Particles)):
        xPositions.append(Particles[i][0])
        yPositions.append(Particles[i][1])
    
    
    xbins=np.arange(0,a+h, h)
    ybins=np.arange(0,a+h, h)
    concentration, xbins, ybins=np.histogram2d(yPositions, xPositions, bins=(xbins, ybins), weights=np.ones_like(xPositions)/(h**2))
    
    return concentration
    

def functionAverage(PreySimulation): 
    
    '''
    Returns the mean-field concentration for each time-step by averaging over all
    simulations
    PreySimulation=list of all Simulations, see Coupling.py
    '''

    timestepsNew=len(PreySimulation[0])
    Average=timestepsNew*[None]
    for t in range(timestepsNew):
        S=[]
        for j in range(len(PreySimulation)):
            
            S.append(Discretization(a, l, PreySimulation[j][t]))
            print(j)
        M=sum(S)/len(S)
        Average[t]=M
        print(len(S),t)
        
    return Average

AllSimulations=[] # collects all simulations
for i in range(30):
    PreySimulation=np.load('./Simulation/Particles'+str(i)+'.npy')
    for j in range(len(PreySimulation)):
        AllSimulations.append(PreySimulation[j])


MeanField=functionAverage(AllSimulations)
np.save('./MeanField', MeanField)

'''Plotting'''

def HybridPlot(Average, Concentration, bd):
	
	'''
	Creates hybrid plots, where the right side corresponds to the FD solution
	and the left side to the mean-field concentration obtained fromt he coupling.
	Average=mean-field concentration
	Concentration=FD solution
	bd=location of the boundary
	'''
	
	listH=[] # list of Hybrid solutions
	for t in range(len(Average)):
		Averaget=Average[t]
		Concentration_t=Concentration[t]
		Hybrid=np.zeros(shape=(l, int(l)))
		for i in range(l):
			for j in range(int(l/2)-bd):
				Hybrid[i,j]=Averaget[i,j]
		
			for j in range(int(l/2)+bd):
				Hybrid[i,j+int(l/2)-bd]=Concentration_t[i,j+int(l/2)-bd]
		listH.append(Hybrid)
	return listH


def functionplot(Data, Max,Time, Name):
	'''
	Creates and saves plots.
	Data=Hybrid or reference solution
	Max=Maximum plotting value
	Time=time-step
	Name=name of figure
	'''

	plt.figure(figsize=(a,a))
	plt.imshow(Data,  interpolation='nearest',cmap='CMRmap', extent=[0, a, 0, a])
	ax=plt.axes()
	cbar=plt.colorbar(fraction=0.045)
	ax.tick_params(labelsize=26)
	cbar.ax.tick_params(labelsize=26)
	plt.xlabel('x', fontsize=26)
	plt.ylabel('y', fontsize=26)
	plt.clim(-Max/20,Max)
	plt.tight_layout()
	plt.savefig('./Plots/Reaction'+str(Time)+str(Name)+'.pdf') 
    
MeanField=np.load('./MeanField.npy')
Reference=np.load('./Simulation/Reference0.npy')
FDSolution=np.load('./FDSolution.npy')
Hybrid=HybridPlot(MeanField, Reference, 0) 

'''Create Plots'''

functionplot(FDSolution[0], 30, 0, 'Reference Initial Condition')


functionplot(Hybrid[4*8], 6, 4, 'Hybrid')
functionplot(Reference[4*8], 6, 4, 'Reference')

functionplot(Hybrid[7*8], 6, 7, 'Hybrid')
functionplot(Reference[7*8], 6, 7, 'Reference')

functionplot(Hybrid[9*8+1], 6, 9, 'Hybrid')
functionplot(Reference[9*8+1], 6, 9, 'Reference')

'''Movie'''  

def functionanimation(Data, Max):
	fig=plt.figure()
	def init():
	 	sns.heatmap(np.zeros((l-1, l-1)), vmin=0, vmax=Max, square=True, cmap='jet', xticklabels=False, yticklabels=False)

	def animate(i):
		data = Data[i]
	
		sns.heatmap(data, annot=False, vmin=0, vmax=Max, cmap='CMRmap', cbar=False, xticklabels=False, yticklabels=False)

	anim = animation.FuncAnimation(fig, animate, init_func=init, frames=int(len(Data)), repeat = False)
	return anim
	
anim=functionanimation(Hybrid, 6)
anim.save('/nfs/datanumerik/bzfkostr/HybridPaper/FirstOrder/Plots/HybridVideo.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
