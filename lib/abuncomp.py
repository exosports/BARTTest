#! /usr/bin/env python

import numpy as np

def compabun(fnamel, fnameh, fnamen, outdir='../results/01BART/'):
    """
    This function takes in 3 spectra (2 with differing abundances, and 1 with 
    the line not present) and computes the factor difference between them.

    Inputs
    ------
    fnamel: string. 'path/to/file' for the spectrum of lesser abundance.
    fnameh: list of strings. 'path/to/file' for the spectra of greater 
            abundances. Note: must be a list even for 1 spectrum.
    fnamen: string. 'path/to/file' for the spectrum without the line. 
    outdir: string. 'path/to/file' for the output text file. Default is for 
            BART.

    Outputs
    -------
    factors: array of the factor differences in abundances.
    05results.txt: Text file of the above factors.

    Example
    -------
    See BARTTest/lib/comparison.py

    Revisions
    ---------
    2017-10-24      mhimes                  Added docstring.
    """
    # Load the spectra
    specl = np.loadtxt(fnamel)
    specn = np.loadtxt(fnamen)
    ninds = np.where(specl[:,1] != np.loadtxt(fnameh[0])[:,1])[0]

    # For all the spectra, store the relation in an array
    factors = np.zeros((len(fnameh), len(ninds)), dtype=float)
    for i in range(len(fnameh)):
        spech = np.loadtxt(fnameh[i])
        # Find where the spectra with lines differ
        inds = np.where(specl[:,1] != spech[:,1])
        # Sum the channels in which the spectra differ
        suml = np.sum(specl[inds])
        sumh = np.sum(spech[inds])
        sumn = np.sum(specn[inds])
        # Compute the factor relation
        factors[i] = ((specn[:,1][inds] - spech[:,1][inds]) / \
                      (specn[:,1][inds] - specl[:,1][inds]))

    # Stats on the factors:
    fmean = np.mean(factors, axis=1)
    fmedn = np.median(factors, axis=1)
    frngl = np.amin(factors, axis=1)
    frngh = np.amax(factors, axis=1)
    fstdv = np.std(factors, axis=1)

    stats = np.stack((fmean, fmedn, frngl, frngh, fstdv), axis=1)

    np.savetxt(outdir+'f05results.txt', stats, fmt='%.6e', \
               header='The following are the factor differences of line ' + \
               'depth between the \n' +                                     \
               'shallowest line (1e-4 abundance) and the others:\n' +       \
               'Mean \t\t Median \t\t Min \t\t Max \t\t Standard deviation')

    return factors


if __name__ == "__main__":
    """
    Produce the plots using the above function
    """
    # Directory with data files
    fdir = '../code-output/01BART/f05abundance/'

    # Spectrum to compare abundances to
    lower =  fdir+'abundance_1e-4_emission_spectrum.dat'
    # List spectra with abundances greater than lower
    upper = [fdir+'abundance_2e-4_emission_spectrum.dat', 
             fdir+'abundance_3e-4_emission_spectrum.dat', 
             fdir+'abundance_4e-4_emission_spectrum.dat', 
             fdir+'abundance_5e-4_emission_spectrum.dat',
             fdir+'abundance_6e-4_emission_spectrum.dat', 
             fdir+'abundance_7e-4_emission_spectrum.dat',
             fdir+'abundance_8e-4_emission_spectrum.dat', 
             fdir+'abundance_9e-4_emission_spectrum.dat',
             fdir+'abundance_1e-3_emission_spectrum.dat']
    norm  =  fdir+'abundance_0_emission_spectrum.dat'

    out = compabun(lower, upper, norm)
    # See the resulting output file (BARTTest/results/01BART/05results.txt) or 
    # the README (BARTTest/tests/05abundance/README) for an explanation of `out`

