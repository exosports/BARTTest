#! /usr/bin/env python

import sys, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs


plt.ion()


fiso   = '../tests/c01hjcleariso/iso_uni.atm'
fnoinv = '../tests/c02hjclearnoinv/noinv_uni.atm'
finv   = '../tests/c03hjclearinv/inv_uni.atm'

iso   = np.loadtxt(fiso,   skiprows=12)[:,1:3]
noinv = np.loadtxt(fnoinv, skiprows=12)[:,1:3]
inv   = np.loadtxt(finv,   skiprows=12)[:,1:3]

fig = plt.figure(tight_layout=True)
grd = gs.GridSpec(1,3)

ax1 = fig.add_subplot(grd[0,0])
ax1.semilogy(iso[:,1], iso[:,0], color='blue', label='isoth')
ax1.set_ylabel('Pressure (bar)')
ax1.set_ylim(np.amax(iso[:,0]), np.amin(iso[:,0]))
ax1.set_xlim(1080, 1180)
ax1.legend(loc='best')

ax2 = fig.add_subplot(grd[0,1])
ax2.semilogy(noinv[:,1], noinv[:,0], color='green', label='non-invert')
ax2.set_xlabel('Temperature (K)')
ax2.set_ylim(np.amax(noinv[:,0]), np.amin(noinv[:,0]))
ax2.legend(loc='best')

ax3 = fig.add_subplot(grd[0,2])
ax3.semilogy(inv[:,1], inv[:,0], color='red', label='invert')
ax3.set_ylim(np.amax(inv[:,0]), np.amin(inv[:,0]))
ax3.legend(loc='upper left')

plt.savefig('../results/01BART/PTs.pdf')
plt.close()







