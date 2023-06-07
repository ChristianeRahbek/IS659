# -*- coding: utf-8 -*-
"""
Created on Fri May  5 14:44:52 2023

@author: Chris
"""

import numpy as np
import matplotlib.pyplot as plt

fileName = 'D:\IS659/finestructure_fit/figures/solidAnglePairData.txt'

data = np.loadtxt(fileName, skiprows=1)

hitAngs = data[:,0]
solidAngs = data[:,1]

merged = [[hitAngs[i], solidAngs[i]] for i in range(0,len(hitAngs))]

sorted(merged, key = lambda angs : merged[0])

sorted(merged, key = lambda k: [k[0],k[1]])

angs = [i[0] for i in merged]
solAngs = [i[1] for i in merged]

fig = plt.figure()
plt.plot(angs, solAngs, ".", markersize=0.1)

fig.savefig('SA_plot.png')

plt.figure()
fig1, ax1 = plt.subplots(figsize = (10,5))

xmin = 0
xmax = 36
hitAngsbins = 5*np.linspace(xmin, xmax, xmax-xmin+1)

solidAngsBins = np.zeros(xmax-xmin+1)


for i in range(xmin, xmax):
    for j in range(len(hitAngs)):
        ang = hitAngs[j]
        if((ang >= hitAngsbins[i]) & (ang < hitAngsbins[i+1])):
            solidAngsBins[i] += solidAngs[j]
    
ax1.bar(hitAngsbins, solidAngsBins, yerr = solidAngsBins*0.01,width = 5, align = 'edge')

mysize = 15

plt.xlabel(r"$\theta_{p_1p_2}$", fontsize = mysize+3)
plt.ylabel(r"$\epsilon_{p_1p_2} / 5\degree$", fontsize = mysize+3)

ax1.ticklabel_format(axis = 'y', style = 'sci', scilimits=(-5,-5))

ax1.tick_params(axis='both', which='both', direction='inout', top = True, right = True, labelsize =mysize)

ax1.yaxis.offsetText.set_fontsize(mysize)

fig1.savefig('SolidAngleEfficiency.png')
print("done")