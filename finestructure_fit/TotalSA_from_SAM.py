# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 09:14:21 2023

@author: Chris
"""

import numpy as np

fileName = 'D:\IS659/finestructure_fit/SAM8He.txt'
#fileName = 'D:\IS659/finestructure_fit/SAM20Na.txt'

U1_SAM = np.loadtxt(fileName, usecols=range(16), skiprows=9, max_rows=16)
U2_SAM = np.loadtxt(fileName, usecols=range(16), skiprows=27, max_rows=16)
U3_SAM = np.loadtxt(fileName, usecols=range(16), skiprows=45, max_rows=16)
U4_SAM = np.loadtxt(fileName, usecols=range(16), skiprows=63, max_rows=16)

U1_SA = U1_SAM.sum()
U2_SA = U2_SAM.sum()
U3_SA = U3_SAM.sum()
U4_SA = U4_SAM.sum()

totalSA = U1_SA + U2_SA + U3_SA + U4_SA

with open("totalSA_8He.txt", "w") as f:
    print('Solid angle for U1: ' + str(U1_SA), file = f)
    print('Solid angle for U2: ' + str(U2_SA), file = f)    
    print('Solid angle for U3: ' + str(U3_SA), file = f)
    print('Solid angle for U4: ' + str(U4_SA), file = f)
    print('Total solid angle : ' + str(totalSA), file = f)
    print('', file = f)
    print('So we cover ' + str(totalSA/(4*np.pi)) + ' of 4pi', file = f)