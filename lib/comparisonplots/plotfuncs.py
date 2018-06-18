import numpy as np
import matplotlib.pyplot as plt

plt.ion()

def readatm(atmfile):
    """
    This function reads an atmospheric file into a Numpy array.

    Inputs
    ------
    atmfile: string. Path/to/file for the atmospheric file used.

    Returns
    -------
    mols   : list of strings. Molecules present in the atmosphere.
    atminfo: array. Shape is (nlayers, nmolecules + 3). For a given layer, 
                    the data stored is radius, pressure, temperature, 
                    mol_1, mol_2, ...
                    Layers are ordered top to bottom. That is, the 0th index 
                    is the top-most layer.

    Example
    -------
    This example assumes that the directory containing Transit is parallel to 
    the current.
    >>> atm = "../transit/transit/examples/demo/HD209458b_demo.atm"
    >>> mols, atminfo, ur, up = atmreader.readatm(atm)

    Revisions
    ---------
    2017-10-02 mhimes@knights.ucf.edu   Initial implementation.
    """
    # Open the atm file and read it
    atmfoo = open(atmfile, 'r')
    lines  = atmfoo.readlines()

    # Store the molecule info
    mols = lines[9].split()

    # Find where the layer data starts
    for i in range(len(lines)):
        line = lines[i].split()
        try:
            if line[0] == '#Radius':
                break
        except:
            continue

    # Trim atm file
    lines = lines[i + 1:]

    # Array to hold atm file info
    atminfo = np.zeros((len(lines), len(mols) + 3), dtype=float)
    # Note that +3 is to hold radius, pressure, and temperature

    # Read in data
    for i in range(len(lines)):
        # Split info for layer `i`
        line = lines[i].split()
        for j in range(3 + len(mols)): # rad, press, temp, mol1, mol2, ...
            # Read in each column
            atminfo[i, j] = float(line[j])

    # Make sure atmosphere is ordered from bottom to top
    if atminfo[0, 0] > atminfo[1, 0]:
        for i in range(atminfo.shape[1]):
            atminfo[:, i] = atminfo[:, i][::-1]

    return mols, atminfo




def makeplots(atmfile, outdirs=['../tests/09comparison/visuals/ptplots/', 
                                '../tests/09comparison/visuals/abunplots/']):
    """
    This function makes a plot of the PT profile and abundance profiles 
    for some atmospheric file.

    Inputs
    ------
    atmifle: string. Path/to/file.ext for the atmospheric file.
    outdirs: list of strings. Paths of where to save the plots, relative to 
             the directory it is being called from.

    Outputs
    -------
    ptplots/atmfile.png  : Plot of the PT profile. File name is the same as 
                           atmfile without the extension.
    abunplots/atmfile.png: Plot of the abundance profiles. File name is the 
                           same as the atmfile without the extension.

    Example
    -------
    [Make an atmospheric file in Transit format named noninverted.atm]
    >>> makeplots('noninverted.atm')

    Revisions
    ---------
    2017-10-16  mhimes@knights.ucf.edu      Initial implementation.
    """
    # Read the atmospheric file
    mols, atminfo = readatm(atmfile)

    # Plot the PT profile and save
    plt.plot(atminfo[:,2], atminfo[:,1], lw=2)
    plt.xlabel('Temperature [K]', fontsize=20)
    plt.ylabel('Pressure [bars]', fontsize=20)
    plt.yscale('log')
    plt.gca().invert_yaxis()
    plt.savefig(outdirs[0] + atmfile.split('/')[-1].split('.')[0] + '_PT.png')
    plt.close()

    # Plot the abundance profiles
    # Don't plot these abundances because they do not contribute opacity
    # Note that HCN is listed: Sharps & Burrows 2007 does not include this as 
    # an opacity source, and RHD uses their opacity table
    ignore = {'H', 'C', 'N', 'O', 'He', 'H2', 'N2', 'HCN', 'CG1', 'CG2', 'CG3'}
    # Set initial values of min and max to extrema
    themin = 1
    themax = 0
    for i in range(len(mols)):
        if mols[i] in ignore:
            continue
        plt.plot(atminfo[:,i+3], atminfo[:,1], label=mols[i], lw=2)
        # Check if this molecule is the minimum or maximum abundance
        if np.amin(atminfo[:,i+3]) < themin:
            themin = np.amin(atminfo[:,i+3])
        if np.amax(atminfo[:,i+3]) > themax:
            themax = np.amax(atminfo[:,i+3])
    plt.xlabel('Molar Mixing Fraction', fontsize=20)
    plt.xscale('log')
    if themin != themax:
        plt.xlim(0.1*themin, 10*themax)
    else:
        plt.xlim(10**-10, 10**-3)
    plt.ylabel('Pressure [bars]', fontsize=20)
    plt.yscale('log')
    plt.gca().invert_yaxis()
    plt.legend(loc='lower left')
    plt.savefig(outdirs[1] + atmfile.split('/')[-1].split('.')[0] + '_abun.png')
    plt.close()


