#!/usr/bin/env python

from openeye import oechem #, oemolbase
import os


class Complex:

    def __init__(self, pdbPath):
        self.pdbPath = pdbPath
        self.residues = {}
        self.ligand = None
        self.pocketResidues = []

        self.correctResNames = ("ALA", "ARG", "ASN", "ASP", "CYS", "GLU",
                                "GLN", "GLY", "HIS", "ILE", "LEU", "LYS",
                                "MET", "PHE", "PRO", "SER", "THR", "TRP",
                                "TYR", "VAL", "HIE", "HSE")

    def makeComplex(self):
        '''
        This function reads the protein, takes each residue, checks if at
        least one of its atoms
        is in the self.consensusResList, if yes it creates a new OEGraphMol
        object representing that residue.
        The OEGraphMol creation is as follows : read each atom,
        read each bond of each atom,
        store all the bonds read, create only when it is a
        new bond, create the atoms in the same way.
        Special care is taken for the N & C terminal of each
        residue, we do not want to select the atom of the
        neibouring residue, a flag helps this careful
        selection.
        '''

        #pdbName = os.path.basename(pdbPath).replace(".pdb", "")
        # Load protein, extract residues
        ifsProt = oechem.oemolistream()
        ifsProt.SetFormat(oechem.OEFormat_PDB)

        if ifsProt.open(self.pdbPath):
            for prot in ifsProt.GetOEGraphMols():

                oechem.OEPerceiveResidues(prot, oechem.OEPreserveResInfo_All)
                hierView = oechem.OEHierView(prot)

                # Looping over the protein's residues
                for mol in hierView.GetResidues():

                    # Check if residue or ligand
                    molName = mol.GetResidueName()
                    if molName in self.correctResNames:

                        # Fix residue number
                        resNumber = str(mol.GetResidueNumber())
                        if len(resNumber) == 1:
                            resNumber = "00" + resNumber
                        elif len(resNumber) == 2:
                            resNumber = "0" + resNumber

                        if molName in ("HIE", "HSE"):
                            molName = "HIS"
                        resTitle = resNumber + "_" + molName

                        # Setup the information to store for this residue
                        if molName == "ASP":
                            residue = self.createMolecule(mol, resTitle, verbose=True)
                        else:
                            residue = self.createMolecule(mol, resTitle)
                        resRings = RingAnalysis(residue).ringsData
                        # Store the information
                        #print resTitle
                        self.residues[resTitle] = [residue, resRings, True]
                    else:
                        # What to store for the ligand
                        ligand = self.createMolecule(mol, molName, verbose=True)
                        ligandRings = RingAnalysis(ligand).ringsData
                        # Store it
                        self.ligand = [ligand, ligandRings]

    def createMolecule(self, mol, molTitle, verbose=False):
        '''
        Get a OE residue, return a OE molecule
        '''

        # Creating a new Res GraphMol
        newMol = oechem.OEGraphMol()
        newMol.SetTitle(molTitle)
        # Create every atom found in the OEHierResidue
        for atom in mol.GetAtoms():
            newMol.NewAtom(atom)

        #---------------------------------------------#
        # Recreate the molecule from atom coordinates #
        #---------------------------------------------#

        # Set the 3D dimension
        #oechem.OESetDimensionFromCoords(newMol)
        # Determine the connectivity between every atom, and assign bond order
        oechem.OEDetermineConnectivity(newMol)
        oechem.OEFindRingAtomsAndBonds(newMol)
        oechem.OEPerceiveBondOrders(newMol)
        # Determine rings

        #oechem.OEAssignAromaticFlags(newMol, oechem.OEAroModelOpenEye)
        # Set hydrogens
        oechem.OEAssignImplicitHydrogens(newMol)
        #oechem.OEAddExplicitHydrogens(newMol)
        # Set charge
        oechem.OEAssignFormalCharges(newMol)

        """
        OEDetermineConnectivity(mol)
        OEFindRingAtomsAndBonds(mol)
        OEPerceiveBondOrders(mol)
        OEAssignImplicitHydrogens(mol)
        OEAssignFormalCharges(mol)
        """


        if verbose:
            print("\nMolecule name: {}\n".format(molTitle))
            for atom in newMol.GetAtoms():
                print(atom, atom.GetFormalCharge())
            ofs = oechem.oemolostream()
            if (ofs.open("test_{}.mol2".format(molTitle)) == 1):
                oechem.OEWriteMolecule(ofs, newMol)

        # Set atom radius
        #oechem.OEAssignBondiVdWRadii(newMol)
        #oechem.OEAssignPartialCharges(newMol,
        #                              oechem.OECharges_None,
        #                              False,
        #                              False)

        return newMol

    def printResidues(self):
        '''
        Print out the residue list
        '''

        sortedKeys = self.residues.keys()
        sortedKeys.sort()

        if len(sortedKeys) == 0:
            print("No residue list loaded")
        else:
            print(','.join(sortedKeys))

    def getResidues(self):
        '''
        Return the list of residue names
        '''
        sortedKeys = sorted(self.residues.keys())

        if len(sortedKeys) == 0:
            print("No residue list loaded")
        else:
            return sortedKeys

    def getResiduesConsensus(self):
        '''
        Return the lst of residues part of the consensus sequence
        '''
        sortedKeys = sorted(self.residues.keys())

        consensusResList = []

        if len(sortedKeys) == 0:
            print("No residue list loaded")
        else:
            for resTitle in sortedKeys:
                residue, resRings, isConsensus = self.residues[resTitle]
                if isConsensus:
                    consensusResList.append(resTitle)

            return consensusResList

    def printLigand(self):
        '''
        Print out the stored ligand
        '''
        if not self.ligand:
            molName = os.path.basename(self.pdbPath)
            print("No ligand loaded in this macromolecule object: " + molName)
        else:
            ligand, ligRings = self.ligand
            print(ligand.GetTitle())

    def getLigand(self):
        '''
        Return the ligand data
        '''
        if not self.ligand:
            molName = os.path.basename(self.pdbPath)
            print("No ligand loaded in this macromolecule object: " + molName)
        else:
            return self.ligand

    """
    def getFprint(self, mode='full'):
        '''
        Extract data from the Complex instance
        '''

        sortedKeys = self.residues.keys()
        sortedKeys.sort()

        resList = []
        fprintList = []

        if len(sortedKeys) == 0:
            print("There is no residue in the molecular complex")
        else:
            for key in sortedKeys:
                residue, resFprint, resRings = self.residues[key]
                # This will return the full fprint sequence, including the
                # residues that show no interaction with the ligand
                if mode == 'full':
                    resList.append(residue.GetTitle())
                    fprintList.append(resFprint)
                # With the optional argument 'stripped' only the residue list
                # that has an fprint that has a interaction with the ligand
                # will be displayed
                elif mode == 'stripped':
                    if "1" in resFprint:
                        resList.append(residue.GetTitle())
                        fprintList.append(resFprint)

        return [resList, fprintList]
    """


class RingAnalysis:
    """
    This class is created to produce ring analysis
    """

    def __init__(self, molecule):
        '''
        The init function calculates information for 5 atom, and 6 atom rings
        (can be extended)
        '''

        # This list is the information accessed from outside this class
        # to get the information
        self.ringsData = []

        # Find aromatic rings
        ringsIdx = []
        ringsIdx += self.determineRings(molecule, 4)
        ringsIdx += self.determineRings(molecule, 5)
        ringsIdx += self.determineRings(molecule, 6)
        ringsIdx += self.determineRings(molecule, 8)

        # Calculate Rings Normal & Center
        for ring in ringsIdx:
            ringData = self.analyseRing(molecule, ring)
            #print ringData
            self.ringsData.append(ringData)

    def determineRings(self, mol, ringSize):
        '''
        This function uses OESbSearch to find all the aromatic rings in a
        given molecule
        It can look for 4 atom, 5 atom, 6 atom and 8 atom aromatic rings
        '''

        # List of rings, each ring being an ordered list of atom idx
        moleculeRings = []

        # Definition of the searches, theses definitions will get sequences of
        # of aromatic atoms, not necessarly rings
        ssAroFour = oechem.OESubSearch("[a;R][a;R][a;R][a;R]")
        ssAroFive = oechem.OESubSearch("[a;R][a;R][a;R][a;R][a;R]")
        ssAroSix = oechem.OESubSearch("[a;R][a;R][a;R][a;R][a;R][a;R]")
        ssAroEight = oechem.OESubSearch("[a;R][a;R][a;R][a;R]\
                                        [a;R][a;R][a;R][a;R]")

        # Adapt search given input
        if ringSize == 4:
            subSearch = ssAroFour
        elif ringSize == 5:
            subSearch = ssAroFive
        elif ringSize == 6:
            subSearch = ssAroSix
        elif ringSize == 8:
            subSearch = ssAroEight

        # match pattern in residue, and cycle through residue atoms
        for patternMatch in subSearch.Match(mol):

            match = patternMatch.GetTargetAtoms()
            ringAtomList = []

            # Store the current match as a list of the atom's Idx
            # and sort the list, so that two lists containing identical
            # atoms lists will be considered as being the same ring
            for atom in match:
                ringAtomList.append(atom.GetIdx())
                ringAtomList.sort()

            # Restore the match iter to the first item,
            # then get the first of the match
            # go to the last item of the match iter
            # then get the last of the match
            match.ToFirst()
            for first in match:
                # going to last ensures the for loop
                # is only used once (no more items to loop over after last)
                match.ToLast()
                # Just store the last OEAtomBase as last
                for last in match:
                    pass

                for bondedToFirst in first.GetAtoms():
                    if bondedToFirst == last:
                        if ringAtomList not in moleculeRings:
                            moleculeRings.append(ringAtomList)
                        # when a cycle is completed, no need to check
                        # the rest of the atoms bonded to first, continue
                        continue

        return moleculeRings

    def analyseRing(self, mol, idxList):
        '''
        This function takes in a list of atom coordinates
        It computes the center of this set of points
        It also computes the normal vector, based on 3 points of the coordList
        (since all those points are in the same plane, being part of
        an aromatic ring)
        It then returns a list containing two lists, those two lists are the
        coordinates of the center, and the coordinates of the normal vector
        (Y-axis) and X-axis and Z-axis

        For each ring, a set of axis is calculated:
            - Y-axis is represented by the normal vector
            - X-axis is the vector between the center of the ring and a
                given atom of the ring
            - Z-axis is the normal to those two previous vectors
        '''

        #----------------#
        # Get Coord list #
        #----------------#

        coordList = []

        for atom in mol.GetAtoms():
            if atom.GetIdx() in idxList:
                coordList.append(mol.GetCoords(atom))

        #--------#
        # Center #
        #--------#

        center = []
        sumX = 0
        sumY = 0
        sumZ = 0

        for count, coord in enumerate(coordList):
            sumX += coord[0]
            sumY += coord[1]
            sumZ += coord[2]

        # Count started at 0, so increment once to have the actual number
        # of a atoms
        count += 1

        center.append(sumX / count)
        center.append(sumY / count)
        center.append(sumZ / count)

        #----------------------#
        # Determine the X-axis #
        #----------------------#

        xAxis = self.defineVector(coordList[0], center)

        #---------------------------------------#
        # Determine Y-axis (normal to the ring) #
        #---------------------------------------#

        radAxis = self.defineVector(coordList[1], center)
        yAxis = self.crossProduct(xAxis, radAxis)

        #------------------#
        # Determine Z-axis #
        #------------------#

        zAxis = self.crossProduct(xAxis, yAxis)

        #---------------#
        # Return result #
        #---------------#

        ringData = []
        ringData.append(center)
        ringData.append(xAxis)
        ringData.append(yAxis)
        ringData.append(zAxis)

        return ringData

    def defineVector(self, coordsA, coordsB):
        '''
        Calculate a vector given two sets of coordinates
        '''

        vector = []

        # Calculating the vectors from the 3 points
        vector.append(coordsA[0] - coordsB[0])
        vector.append(coordsA[1] - coordsB[1])
        vector.append(coordsA[2] - coordsB[2])

        return vector

    def crossProduct(self, vecA, vecB):
        '''
        Determine the cross product between two vectors
        (the perpendicular to the plane formed by those two vectors)
        Formula:
        A & B : [   Ay*Bz - Az*By   ,   Az*Bx - Ax*Bz   ,   Ax*By - Ay*Bx   ]
        '''

        crossProd = []
        crossProd.append(vecA[1] * vecB[2] - vecA[2] * vecB[1])
        crossProd.append(vecA[2] * vecB[0] - vecA[0] * vecB[2])
        crossProd.append(vecA[0] * vecB[1] - vecA[1] * vecB[0])

        return crossProd
