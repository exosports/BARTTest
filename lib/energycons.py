#! /usr/bin/env python

import sys
import numpy as np
#import scipy.interpolate as si
sys.path.append("../../BART/modules/transit/scripts/")
import readtransit as rt

def energycons(spectra, outfile):
    """
    Integrates each spectrum, saves those values as well as if they all 
    match within a specified tolerance.

    Inputs
    ------
    spectra: list, strings. Path/to/file for each spectrum to be integrated.
    outfile: string.        Path/to/file to save the results.
    """
    # Array to hold integrated values
    integspec = np.zeros(len(spectra))

    # Load each spectrum, integrate it
    for i in range(len(spectra)):
        wnums, flux  = rt.readspectrum(spectra[i])
        integspec[i] = np.trapz(flux, x=wnums)
        # Cubic spline returns roughly same value
        #integspec[i] = si.InterpolatedUnivariateSpline(wnums, flux, 
        #                                 k=3).integral(wnums[0], wnums[-1])

    # Check if they match
    diff = integspec / integspec[0] - 1

    # Save out result
    hdr = "# Units are erg/s/cm2"
    ftr = "# % differences compared to the first spectrum: \n" + \
          "#   " + str(100*diff) + " %"
    np.savetxt(outfile, integspec, header=hdr, footer=ftr)


if __name__ == "__main__":
    spectra = ['../code-output/01BART/f09energycons/' + \
                'energycons_1_emission_spectrum.dat',
               '../code-output/01BART/f09energycons/' + \
                'energycons_2_emission_spectrum.dat',
               '../code-output/01BART/f09energycons/' + \
                'energycons_5_emission_spectrum.dat',
               '../code-output/01BART/f09energycons/' + \
                'energycons_10_emission_spectrum.dat'   ]
    energycons(spectra, '../results/01BART/f09energycons.txt')






