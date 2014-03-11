#!/usr/bin/env python

__description__ = \
    """
    PDB_renumSplitWaterChain.py

    Takes a .pdb water chain produced by VMD and splits it into as many .pdb
    files as needed not to exceed the max 'residue number'
    """
__author__ = "Thomas Coudrat"
__date__ = "131019"


import sys
import os


def main():
    """
    Run this script
    """

    # Get the input
    inputPath = sys.argv[1]
    inputFile = open(inputPath, 'r')
    inputLines = inputFile.readlines()
    inputFile.close()

    # Get the split of lines, and header and tail
    lineGroups, headLine, tailLine = get_renumbered_lines(inputLines)

    inputName = os.path.basename(inputPath)
    for i, newLines in enumerate(lineGroups):

        # Use the input file name with a suffix to name the output file names
        outName = inputName.replace(".pdb", "_" + str(i + 1) + ".pdb")
        outFile = open(outName, "w")

        # Write the lines, including header and tail
        outFile.write(headLine)
        for line in newLines:
            outFile.write(line)
        outFile.write(tailLine)


def get_renumbered_lines(inputLines):
    """
    Read the lines of the input file
    Modify the ATOM lines to renumber the atom and residue number
    Return those renumbered lines
    Also return the head and tail lines of the input file
    """

    # Store lines to return here
    headLine = inputLines[0]
    tailLine = inputLines[-1]
    # In this store list of lists of lines
    lineGroups = []

    # New numbers for residues and atoms
    newResNum = 0
    newAtomNum = 0
    # This just keeps track of the previous residue number when looping
    prevResNum = 0

    for line in inputLines:
        l = line.split()
        if l[0] == 'ATOM' or l[0] == 'HETATM':
            lineType = line[0:6].strip()
            #atmNum = line[6:11].strip()
            atmType = line[12:17].strip()
            resName = line[17:21].strip()
            chain = line[21].strip()
            currResNum = int(line[22:26])
            x = line[30:38].strip()
            y = line[38:46].strip()
            z = line[46:54].strip()
            occup = line[54:60].strip()
            tempFac = line[60:66].strip()
            elem = line[72:75].strip()

            # Update the residue number to +1 compared to previous res number
            if currResNum != prevResNum:
                prevResNum = currResNum
                newResNum += 1

            # Update atom number
            newAtomNum += 1

            # When the new residue number goes over the max, change it back to
            # the number 0. Also reset newAtomNum
            if newResNum > 9999:
                newResNum = 1
                newAtomNum = 1

            # Create a new lineGroup instance every time the newResNum = 1, and
            # newAtomNum = 1 which is at first loop turn, and everytime
            # the newResNum reaches more than 9999
            if newResNum == 1 and newAtomNum == 1:
                lineGroups.append([])
                currentLines = lineGroups[-1]

            # pdb format line
            newline = '%-5s%6s  %-3s%5s%s%4s%12s%8s%8s%6s%6s%9s' \
                % (lineType, newAtomNum, atmType, resName, chain, newResNum,
                   x, y, z, occup, tempFac, elem)
            #print newline
            currentLines.append(newline + "\n")

    return lineGroups, headLine, tailLine


if __name__ == "__main__":
    main()
