#! /usr/bin/env python
"""
Script used to create Barstow et al. CO-only atmospheric models
This script was used in two parts
The first part is the commented out block of code below the first set of imports
That code creates a plaintext file of the atmospheric model values, ordered 
pressure, temperature, gas abun 1, gas abun 2, ...
After that, the TEA header was added to the files
That code was then commented out, and the code below that was executed to 
calculate the radius at each pressure, temperature point
"""

import numpy as np
import sys
import scipy.interpolate as sci
import matplotlib.pyplot as plt
'''
plog = np.linspace(-5,1,num=100)
PT = np.zeros((100,5))
PT[:,0] = 10**plog
PT[:,1] = 1000
PT[:,4] = 1e-4
PT[:,2] = 0.85 / (0.85+0.15) * (1 - 1e-4)
PT[:,3] = 0.15 / (0.85+0.15) * (1 - 1e-4)
np.savetxt('CO_1e-4_1000K2.atm', PT, fmt='%.09e')
PT[:,1] = 1500
np.savetxt('CO_1e-4_1500K2.atm', PT, fmt='%.09e')
PT[:,4] = 1e-5
PT[:,2] = 0.85 / (0.85+0.15) * (1 - 1e-5)
PT[:,3] = 0.15 / (0.85+0.15) * (1 - 1e-5)
np.savetxt('CO_1e-5_1500K2.atm', PT, fmt='%.09e')
'''

sys.path.append('../../../BART/code/')
import makeatm as ma

outspec = 'H2_ref He_ref CO_g'

abun_file = '../../../BART/inputs/abundances_Asplund2009.txt'
refpress = 10.0 #bars, Tinetti et al. (2012), ref'd by Barstow et al. (2020)
tep_name = '../00inputs/BarstowEtal_CO.tep'
ma.makeRadius(outspec, 'CO_1e-4_1000K2.atm', abun_file, tep_name, refpress)
ma.makeRadius(outspec, 'CO_1e-4_1500K2.atm', abun_file, tep_name, refpress)
ma.makeRadius(outspec, 'CO_1e-5_1500K2.atm', abun_file, tep_name, refpress)

