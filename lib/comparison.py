#! /usr/bin/env python

"""
This file computes the factor differences in line depth for a number of 
different spectra. For more details, see the docstring of compabun as this 
file serves more as an example of that.
"""

from abuncomp import compabun

def main():
    # Directory with data files
    fdir = '../code-output/01BART/'

    # Spectrum to compare abundances to
    lower =  fdir+'05abundance_1e-4_emission_spectrum.dat'
    # List spectra with abundances greater than lower
    upper = [fdir+'05abundance_2e-4_emission_spectrum.dat', 
             fdir+'05abundance_3e-4_emission_spectrum.dat', 
             fdir+'05abundance_4e-4_emission_spectrum.dat', 
             fdir+'05abundance_5e-4_emission_spectrum.dat',
             fdir+'05abundance_6e-4_emission_spectrum.dat', 
             fdir+'05abundance_7e-4_emission_spectrum.dat',
             fdir+'05abundance_8e-4_emission_spectrum.dat', 
             fdir+'05abundance_9e-4_emission_spectrum.dat',
             fdir+'05abundance_1e-3_emission_spectrum.dat']
    norm  =  fdir+'05abundance_linemoved_emission_spectrum.dat'

    out = compabun(lower, upper, norm)
    # See the resulting output file (BARTTest/results/01BART/05results.txt) or 
    # the README (BARTTest/tests/05abundance/README) for an explanation of `out`

if __name__ == "__main__":
    main()

