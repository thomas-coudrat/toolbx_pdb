#!/usr/bin/python

# Takes a directory containing .pdb files as argument, and computes the average
# b-factor of all heavy atoms with the ATOM prefixes in these .pdb files
#
# https://github.com/thomas-coudrat/toolbx_pdb
# Thomas Coudrat <thomas.coudrat@gmail.com>

import argparse
import sys
from glob import glob

def main():
    """
    Run the bfac script
    """

    # Get the first argument of the executed script
    pdbDir = parsing()

    # Print the Bfactors from the X-ray structure files in this directory
    print_bfactors(pdbDir)


def print_bfactors(pdbDir):
    """
    Read each PDB file, collect B-factor values and atom count.
    Calculate average B-factor.
    """

    print("\nAnalysing B-factor from X-ray structures in directory: {}\n".format(pdbDir))

    # Make a list of all .pdb files in the directory given as argument
    pdbFilePaths = glob(pdbDir + "/*.pdb")

    # Loop over all those .pdb files to extract the information
    for pdbPath in pdbFilePaths:

        # Read the pdb
        pdbFile = open(pdbPath, "r")
        pdbLines = pdbFile.readlines()
        pdbFile.close()

        # Get the b-factor list for all atoms in that list
        bFactors = []
        for line in pdbLines:
            ll = line.split()
            # if the line represents an ATOM, and if this atom is not a hydrogen
            if ll[0] == 'ATOM' and ll[2][0] != 'H':
                # b-factors are written in the pdb files between columns 60 and 67
                bFac = line[60:67].strip()
                bFac = float(bFac)
                bFactors.append(bFac)

        # Loop over the b-factors, and get the total and average
        bFacTotal = 0
        numAtoms = len(bFactors)
        for line in bFactors:
            bFacTotal = bFacTotal + line

        print("B-factor on:\t {}".format(pdbPath))
        print("total bFac:\t {}".format(bFacTotal))
        print("atom count:\t {}".format(numAtoms))
        print("bFac avg:\t {}\n".format(bFacTotal / numAtoms))

def parsing():
    """
    Parse arguments
    """

    descr = "Display the B-factor of X-ray structure PDB files"
    descr_pdbDir = "Directory containing X-ray structure files"

    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument("pdbDir", help=descr_pdbDir)
    args = parser.parse_args()
    pdbDir = args.pdbDir

    return pdbDir

if __name__ == "__main__":
    main()
