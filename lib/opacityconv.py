import numpy as np
import struct

def bin2np(foo, output='opacity'):
    """
    This function takes an opacity file (binary) produced by transit and 
    converts it to Numpy arrays.

    Inputs
    -------
    foo   : String. Transit opacity file.
    output: String. path/to/filename for saved Numpy arrays. Do not include 
                    extension, but DO include a file name. 
                    Default file name is opacity(.npz)

    Outputs
    -------
    header : Numpy array. Nmol, Ntemp, Nlayer, Nwave values.
    molID  : Numpy array. Molecule isotopologue.
    temp   : Numpy array. Array of temperatures for opacity table.
    press  : Numpy array. Array of pressures for opacity table.
    wns    : Numpy array. Array of wavenumbers for opacity table.
    optable: Numpy array. Opacity table.

    Example
    -------
    >>> header, molID, temp, press, wns, optable = bin2np('./opacity.opt')
    For the creation of `opacity.opt`, run the broadening test for BART.

    Revisions
    ---------
    2017-10-24      mhimes                  Added example and revisions to doc
    """
    # Load all binary data
    with open(foo, 'rb') as foop:
        filecontent = foop.read()

    # Read header values: Nmol, Ntemp, Nlayer, Nwave
    hlen = 4*8 #4 values, 8 bytes each
    header = struct.unpack("llll", filecontent[:hlen])
    header = np.asarray(header)
    Nmol, Ntemp, Nlayer, Nwave = header

    # Read in mol ID array
    mlen = hlen + Nmol*4 #Nmol values, 4 bytes each
    molID = struct.unpack("i" * Nmol, filecontent[hlen:mlen])
    molID = np.asarray(molID)

    # Read in temp array
    tlen = mlen + Ntemp*8 #Ntemp values, 8 bytes each
    temp = struct.unpack("d" * Ntemp, filecontent[mlen:tlen])
    temp = np.asarray(temp)

    # Read in pressure array
    plen = tlen + Nlayer*8 #Nlayer values, 8 bytes each
    press = struct.unpack("d" * Nlayer, filecontent[tlen:plen])
    press = np.asarray(press)
    # All values are 1e6 larger than actual
    press = press * 1e-6

    # Read in wavenumber array
    wlen = plen + Nwave*8 #Nwave values, 8 bytes each
    wns = struct.unpack("d" * Nwave, filecontent[plen:wlen])
    wns = np.asarray(wns)

    # Read in opacity table, then reshape
    optable = struct.unpack("d"*Nwave*Nmol*Ntemp*Nlayer, filecontent[wlen:])
    optable = np.asarray(optable).reshape(Nlayer, Ntemp, Nmol, Nwave)

    # Save arrays to Numpy file
    np.savez(output+'.npz', header=header, molID=molID, temp=temp,            \
             press=press, wns=wns, optable=optable)

    # Return arrays
    return (header, molID, temp, press, wns, optable)

