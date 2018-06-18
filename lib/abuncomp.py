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

    # For all the spectra, store the relation in an array
    factors = np.zeros(len(fnameh), dtype=float)
    for i in range(len(fnameh)):
        spech = np.loadtxt(fnameh[i])
        # Find where the spectra with lines differ
        inds = np.where(specl[:,1] != spech[:,1])
        # Sum the channels in which the spectra differ
        suml = np.sum(specl[inds])
        sumh = np.sum(spech[inds])
        sumn = np.sum(specn[inds])
        # Compute the factor relation
        factors[i] = (sumn - sumh)/(sumn - suml)

    np.savetxt(outdir+'05results.txt', factors, fmt='%.6e', \
               header='The following are the factor differences of line ' + \
               'depth between the \n' +                                     \
               'shallowest line (1e-4 abundance) and the others:')

    return factors
