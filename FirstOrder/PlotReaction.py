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
import animation

'''Parameters'''

l=100 # nuber of cells
a=12 # length of the domain
D=0.5 # Diffusion constant
L=a/2

dx_hist=a/l
deltat=0.01 # Time-step size

'''Discretization'''

def functiondiscretization(a, discret, Particles): # discretizes 
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
    

def functionAverage(PreySimulation): # averages over all simulations

    timestepsNew=len(PreySimulation[0])
    Average=timestepsNew*[None]
    for t in range(timestepsNew):
        S=[]
        for j in range(len(PreySimulation)):
            
            S.append(functiondiscretization(a, l, PreySimulation[j][t]))
  
        M=sum(S)/len(S)
        Average[t]=M
        
    return Average

'''Plotting'''

def functionHybridPlot(Average, ConcentrationData, bd): # put continuous and particle plots together
	listH=[]
	for t in range(len(Average)):
		Averaget=Average[t]
		ConcentrationDatat=ConcentrationData[t]
		H=np.zeros(shape=(l, int(l)))
		for i in range(l):
			for j in range(int(l/2)-bd):
				H[i,j]=Averaget[i,j]
		
			for j in range(int(l/2)+bd):
				H[i,j+int(l/2)-bd]=ConcentrationDatat[i,j+int(l/2)-bd]
		listH.append(H)
	return listH


def functionplot(Data, Max,Time, Name):

	plt.figure(figsize=(12,12))
	plt.imshow(Data,  interpolation='nearest',cmap='CMRmap', extent=[0, 12, 0, 12])
	ax=plt.axes()
	cbar=plt.colorbar(fraction=0.045)
	ax.tick_params(labelsize=26)
	cbar.ax.tick_params(labelsize=26)
	plt.xlabel('x', fontsize=26)
	plt.ylabel('y', fontsize=26)
	plt.clim(-Max/20,Max)
	plt.tight_layout()
	plt.savefig('/data/numerik/bzfkostr/HybridPaper/FirstOrder/Plots/Reaction'+str(Time)+str(Name)+'.pdf') 
    
Average=np.load('/data/numerik/bzfkostr/HybridPaper/FirstOrder/Data/ReactionAverage.npy')
Reference=np.load('/data/numerik/bzfkostr/HybridPaper/FirstOrder/Data/ReactionReference0.npy')
listC=np.load('/data/numerik/bzfkostr/HybridPaper/FirstOrder/Data/1501F.npy')
Hybrid=functionHybridPlot(Average, Reference, 0) 

'''Create Plots'''

functionplot(listC[0], 30, 0, 'Reference Initial Condition')


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
