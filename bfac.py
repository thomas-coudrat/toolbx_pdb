#!/usr/bin/python

#
# Takes a directory containing .pdb files as argument, and computes the average
# b-factor of all heavy atoms with the ATOM prefixes in these .pdb files
#
# Thomas Coudrat, July 2013
#

import sys
from glob import glob

# Get the first argument of the executed script
pdbDirPath = sys.argv[1:][0]

# Make a list of all .pdb files in the directory given as argument
pdbFilePaths = glob(pdbDirPath + "/*.pdb")

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

    print "B-factor on:\t", pdbPath
    print "total bFac:\t", bFacTotal
    print "atom count:\t", numAtoms
    print "bFac avg:\t", bFacTotal / numAtoms
    print
