#! /usr/bin/env python
"""
This script created some of the files in this directory
Originally, this test was only conducted for wavelengths > 1 micron, 
but later it was changed to include the full range of wavelengths considered 
by Barstow et al.
As a result, the files produced by this script were renamed to have 
an "_over1um" suffix, to correspond to the "_under1um" suffix adopted for 
the files created after the test was changed.
This script has not been changed in light of that, and is meant as a reference 
if someone were interested in how the files were created

Note: The input to this script comes from the Barstow et al. online compendium
      The files are named, <code>_hd189733b-cloudfree_jwst_60.dat
These files were previously split to be above/under 1um, to correspond to the 
chosen wavelength region for the test.  This split in the input files has been 
left for historical reasons; if one were to repeat this analysis now, there 
would be no reason for that split (see the s05 test, for example).
"""

import sys
import numpy as np


foo     = sys.argv[1]
midbins = np.loadtxt(foo)
code    = foo.split('_')[0]

# Save out mdibins, depths, uncs
np.savetxt(code+'_midbins.txt', midbins[:,0])
np.savetxt(code+'_depths.txt', midbins[:,1])
np.savetxt(code+'_uncs.txt', midbins[:,2])

midbins = midbins[:,0]


eps = 1e-9

lftbins = np.zeros(len(midbins))
rgtbins = np.zeros(len(midbins))
lftbins[ 0 ] = midbins[ 0 ] - (midbins[ 1] - midbins[ 0 ])/2.
lftbins[ 1:] = midbins[ 1:] - (midbins[1:] - midbins[:-1])/2.
rgtbins[-1 ] = midbins[-1 ] + (midbins[-1] - midbins[-2 ])/2.
rgtbins[:-1] = midbins[:-1] + (midbins[1:] - midbins[:-1])/2.

fnames = []
for i in range(len(midbins)):
    fname = '../00inputs/filters/'+code+'_'+str(i+1).zfill(3)+'.dat'
    with open(fname, 'w') as foop:
        foop.write(str(lftbins[i] - 2*eps)         + ' 0.000\n')
        foop.write(str(lftbins[i] -   eps)         + ' 0.000\n')
        foop.write(str(lftbins[i]        )         + ' 1.000\n')
        foop.write(str((lftbins[i]+rgtbins[i])/2.) + ' 1.000\n')
        foop.write(str(rgtbins[i]        )         + ' 1.000\n')
        foop.write(str(rgtbins[i] +   eps)         + ' 0.000\n')
        foop.write(str(rgtbins[i] + 2*eps)         + ' 0.000\n')
    fnames.append(fname)

with open(code+'_filters.txt', 'w') as foop:
    for i in range(len(fnames)):
        foop.write(fnames[i] + '\n')

