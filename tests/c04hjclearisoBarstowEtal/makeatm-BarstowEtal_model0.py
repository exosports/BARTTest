#! /usr/bin/env python
"""
Script used to create Barstow et al. example hot Jupiter atmospheric model
This script was used in two parts
The first part is the commented out block of code below the first set of imports
That code creates a plaintext file of the atmospheric model values, ordered 
pressure, temperature, gas abun 1, gas abun 2, ...
After that, the TEA header was added to the file
That code was then commented out, and the code below that was executed to 
calculate the radius at each pressure, temperature point
"""

import numpy as np
import sys
import scipy.interpolate as sci
import matplotlib.pyplot as plt
'''
plog = np.linspace(-5,1,num=100)
PT = np.zeros((100,6))
PT[:,0] = 10**plog
PT[:,1] = 1500
PT[:,4] = 0.0003
PT[:,5] = 0.00035
PT[:,2] = 0.85 / (0.85+0.15) * (1 - 0.0003 - 0.00035)
PT[:,3] = 0.15 / (0.85+0.15) * (1 - 0.0003 - 0.00035)
np.savetxt('model0.atm', PT, fmt='%.09e')
'''

sys.path.append('../../../BART/code/')
import makeatm as ma

outspec = 'H2_ref He_ref H2O_g CO_g'

abun_file = '../../../BART/inputs/abundances_Asplund2009.txt'
refpress = 10.0 #bars, Tinetti et al. (2012), ref'd by Barstow et al. (2020)
tep_name = '../00inputs/BarstowEtal_model0.tep'
ma.makeRadius(outspec, 'model0.atm', abun_file, tep_name, refpress)

