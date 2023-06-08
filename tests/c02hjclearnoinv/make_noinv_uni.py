#! /usr/bin/env python

import sys
import numpy as np
import astropy.constants as ac

sys.path.append('../../../BART/code/')
import PT


press = np.logspace(-8, 2, 100)
kappa  = -1.0
gamma1 = -1.3
gamma2 = -0.1
alpha  =  0.6
beta   =  0.8
R_star =  0.756 * ac.R_sun.value
T_star = 5000
T_int  = 100
sma    = 0.03099 * ac.au.value
grav   = 2182.73
T_int_type = 'const'

temp = PT.PT_line(press, kappa,  gamma1, gamma2, alpha, beta, R_star, 
                  T_star, T_int, sma, grav, T_int_type)

out = np.zeros((100, 9))
out[:,0] = press
out[:,1] = temp
out[:,2] = 0.1499751
out[:,3] = 0.8498589
out[:,4] = 1e-4
out[:,5] = 1e-5
out[:,6] = 5e-6
out[:,7] = 5e-5
out[:,8] = 1e-6
np.savetxt('noinv_uni.atm', out, fmt='%.8e', comments='', header='''# This is a final TEA output file with calculated abundances (mixing fractions) for all listed species.
# Units: pressure (bar), temperature (K), abundance (unitless).

#Values units:
ur 1e5
up 1e6
q number

#SPECIES
He H2 CO CO2 CH4 H2O NH3

#TEADATA
#Radius    Pressure   Temp    He         H2         CO         CO2        CH4        H2O        NH3           ''')

import makeatm as ma

outspec = 'He_ref H2_ref CO_g CO2_g CH4_g H2O_g NH3_g'

abun_file = '../../../BART/inputs/abundances_Asplund2009.txt'
refpress = 0.1 #bars
tep_name = '../00inputs/HD189733b.tep'
ma.makeRadius(outspec, 'noinv_uni.atm', abun_file, tep_name, refpress)

print("Now update noinv_uni.atm's header to match that in make_noinv_uni.py")
