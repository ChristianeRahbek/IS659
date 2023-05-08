# -*- coding: utf-8 -*-
"""
Created on Fri May  5 14:44:52 2023

@author: Chris
"""

import numpy as np
import matplotlib.pyplot as plt

fileName = 'D:\IS659/finestructure_fit/figures/solidAnglePairData.txt'

data = np.loadtxt(fileName, skiprows=1)

plt.plot(data[:,0], data[:,1], ".")