import numpy as np
from scipy.special import wofz

"""
From http://scipython.com/book/chapter-8-scipy/examples/the-voigt-profile/
"""

def G(x, alpha, shift=0):
    """ Return Gaussian line shape at x with HWHM alpha """
    return np.sqrt(np.log(2) / np.pi) / alpha\
                             * np.exp(-((x - shift) / alpha)**2 * np.log(2))

def L(x, gamma, shift=0):
    """ Return Lorentzian line shape at x with HWHM gamma """
    return gamma / np.pi / ((x - shift)**2 + gamma**2)

def V(x, alpha, gamma, shift=0):
    """
    Return the Voigt line shape at `x` + `shift` with Lorentzian component 
    HWHM gamma and Gaussian component HWHM alpha.

    """
    sigma = alpha / np.sqrt(2 * np.log(2))

    return np.real(wofz(((x - shift) + 1j*gamma)/sigma/np.sqrt(2))) / sigma\
                                                           /np.sqrt(2*np.pi)


