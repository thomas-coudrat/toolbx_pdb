#!/usr/bin/env python

import argparse
import os
import sys

# Argument parsing
parser = argparse.ArgumentParser(description="Reads a .pdb file and creates as\
    many new .pdb files as needed to overcome the max number of residue\
    reached in a single chain")
parser.add_argument("inputPDB", help="Input .pdb file")
args = parser.parse_args()
pdbPath = args.inputPDB

# Read the lines of the input file
pdbFile = open(pdbPath, 'r')
pdbLines = pdbFile.readlines()
pdbFile.close()

# Extract head and tail to be added to each new .pdb file
pdbHead = pdbLines[0]
pdbTail = pdbLines[-1]

# This dictionary will hold the lines of each new .pdb file, the key of the
# dictionary is the 'WTx' number
waterDict = {}

# Go through the lines of ATOM in the input pdb file
for line in pdbLines[1:-1]:
    # Grab the WTx number which defines when the numbering is restarted as
    # the residue number has reached its maximum
    l = line.split()
    waterGroup = l[-1]

    # Add each line to its corresponding list, based on the key, create a new
    # list when the key was not present yet
    if waterGroup not in waterDict:
        waterDict[waterGroup] = [line]
    else:
        waterDict[waterGroup].append(line)

if len(waterDict) == 1:
    print "This file does not exceed the max residue numbers per chain"
    sys.exit()

# Loop over the dictionary created
print "The following files will get created:"
for waterGroup in waterDict:
    waterNumber = waterGroup.replace("WT", "")
    newPdbPath = pdbPath.replace(".pdb", "_" + waterNumber + ".pdb")
    atomCount = len(waterDict[waterGroup])
    waterCount = atomCount / 3
    print os.path.basename(newPdbPath),\
        "\t atom count:", atomCount,\
        "\t water count:", waterCount

# Prompt for user input
answer = raw_input("Do you want to proceed with file creation? (y/n) ").strip()

if answer == "y":
    # Loop over the dictionary created
    for waterKey in waterDict:
        # Create a new waterPdb file
        waterNumber = waterKey.replace("WT", "")
        waterPath = pdbPath.replace(".pdb", "_" + waterNumber + ".pdb")
        waterFile = open(waterPath, "w")

        # Loop over the stored lines and write that new file
        waterLines = waterDict[waterKey]
        waterFile.write(pdbHead)
        for waterLine in waterLines:
            waterFile.write(waterLine)
        waterFile.write(pdbTail)
        waterFile.close()

elif answer == "n":
    print "Ok, goodbye!"
else:
    print "Please use type either y or n"
