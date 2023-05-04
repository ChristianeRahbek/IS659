# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 09:14:21 2023

@author: Chris
"""

import numpy as np

fileName = 'D:\IS659/finestructure_fit/SAM8He.txt'
#fileName = 'D:\IS659/finestructure_fit/SAM20Na.txt'

fileName_pUS = 'D:\IS659/finestructure_fit/SAM8He_pUS.txt'
fileName_mUS = 'D:\IS659/finestructure_fit/SAM8He_mUS.txt'

U1_SAM = np.loadtxt(fileName, usecols=range(16), skiprows=9, max_rows=16)
U2_SAM = np.loadtxt(fileName, usecols=range(16), skiprows=27, max_rows=16)
U3_SAM = np.loadtxt(fileName, usecols=range(16), skiprows=45, max_rows=16)
U4_SAM = np.loadtxt(fileName, usecols=range(16), skiprows=63, max_rows=16)

U1pUS_SAM = np.loadtxt(fileName_pUS, usecols=range(16), skiprows=9, max_rows=16)
U2pUS_SAM = np.loadtxt(fileName_pUS, usecols=range(16), skiprows=27, max_rows=16)
U3pUS_SAM = np.loadtxt(fileName_pUS, usecols=range(16), skiprows=45, max_rows=16)
U4pUS_SAM = np.loadtxt(fileName_pUS, usecols=range(16), skiprows=63, max_rows=16)

U1mUS_SAM = np.loadtxt(fileName_mUS, usecols=range(16), skiprows=9, max_rows=16)
U2mUS_SAM = np.loadtxt(fileName_mUS, usecols=range(16), skiprows=27, max_rows=16)
U3mUS_SAM = np.loadtxt(fileName_mUS, usecols=range(16), skiprows=45, max_rows=16)
U4mUS_SAM = np.loadtxt(fileName_mUS, usecols=range(16), skiprows=63, max_rows=16)


U1_SA = U1_SAM.sum()
U2_SA = U2_SAM.sum()
U3_SA = U3_SAM.sum()
U4_SA = U4_SAM.sum()

U1pUS_SA = U1pUS_SAM.sum()
U2pUS_SA = U2pUS_SAM.sum()
U3pUS_SA = U3pUS_SAM.sum()
U4pUS_SA = U4pUS_SAM.sum()

U1mUS_SA = U1mUS_SAM.sum()
U2mUS_SA = U2mUS_SAM.sum()
U3mUS_SA = U3mUS_SAM.sum()
U4mUS_SA = U4mUS_SAM.sum()

totalSA = U1_SA + U2_SA + U3_SA + U4_SA
totalSAp = U1pUS_SA + U2pUS_SA + U3pUS_SA + U4pUS_SA
totalSAm = U1mUS_SA + U2mUS_SA + U3mUS_SA + U4mUS_SA

totalSApUS = np.abs(totalSA-totalSAp)
totalSAmUS = np.abs(totalSA-totalSAm)


with open("totalSA_8He.txt", "w") as f:
    print('Solid angle for U1: ' + str(U1_SA), file = f)
    print('Solid angle for U2: ' + str(U2_SA), file = f)    
    print('Solid angle for U3: ' + str(U3_SA), file = f)
    print('Solid angle for U4: ' + str(U4_SA), file = f)
    print('Total solid angle : ' + str(totalSA), file = f)
    print('Plus uncertainty of total solid angle: ' + str(totalSApUS), file = f)
    print('Minus uncertainty of total solid angle: ' + str(totalSAmUS), file = f)
    print('', file = f)
    print('So we cover ' + str(totalSA/(4*np.pi)) + '+- ' + str(totalSAmUS/(4*np.pi)) +' of 4pi', file = f)
    print('', file = f)
    print('The error is ' + str(totalSAmUS/totalSA*100) + '% of the calculated TSA.', file = f)