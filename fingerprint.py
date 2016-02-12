#!/usr/bin/env python


import math
from openeye import oechem  # , oegrid


class Fingerprint:

    def __init__(self, molComplex, fprintDef=False):
        self.molComplex = molComplex
        self.fprint = []
        self.fprintConsensus = []
        if fprintDef == None:
            self.fprintDef = "11111111111"
        else:
            self.fprintDef = fprintDef

        #---------------#
        # OESubSearches #
        #---------------#

        # Definition of a Hydrophobe atom
        # Any atom (Carbon, or Sulfur, or Fluor, or Chloride, or Bromine,
        # or Iodine)
        self.ssHydrophobe = oechem.OESubSearch('[C,S,F,Cl,Br,I]')

        # Definition of a H-Bond donor
        # All atoms (Oxygen, or Nitrogen, or Sulfur, or Fluor)
        # connected to (a Hydrogen)
        self.ssDonor = oechem.OESubSearch('[O,N,S,F][H]')

        # Definition of a H-Bond acceptor
        # All atoms (Oxygen, or Nitrogen, or any negatively charged atom)
        # AND (not a positively charged atom)
        self.ssAcceptor = oechem.OESubSearch('[O,N,*-;!+]')

        # Definition of a weak H-bond Acceptor
        # Any 2 Aromatic atoms connected by aromatic bond,
        # or 2 Aliphatic atoms connected by double bond,
        # or 2 Aliphatic atoms connected by triple bond,
        # self.ssWkAcceptor = OESubSearch("[a:a, A=A, A#A]")
        self.ssWkAcceptor = oechem.OESubSearch('[A,a]=,#,:[A,a]')

        # Definition of a weak H-bond Donor
        # Any (Aromatic carbon atom, or Aliphatic carbon with 3 total bonds
        # or Aliphatic carbon with 2 total bonds) connected to (a Hydrogen)
        self.ssWkDonor = oechem.OESubSearch('[c,CX3,CX2][H]')

        # Definition of a Cation :: Any negatively charged atom
        self.ssCation = oechem.OESubSearch('[*+]')

        # Definition of a Anion :: Any positively charged atom
        self.ssAnion = oechem.OESubSearch('[*-]')

        # Definition of a metal atom
        self.ssMetal = oechem.OESubSearch('[Ca,Cd,Co,Cu,Fe,Mg,Mn,Ni,Zn]')

    def vectorAngle(self, vecAB, vecCD):
        '''
        This function calculates the angle between two vectors.
        The input is 2 vectors, each represented by its x,y,z coordinates.
        The returned value is the angle in radians between vector AB, and
        vector CD.
        '''

        # Get dotProduct of the vectors, and magnitude of each vector
        dotProduct = \
            vecAB[0] * vecCD[0] + \
            vecAB[1] * vecCD[1] + \
            vecAB[2] * vecCD[2]
        magnitudeAB = math.sqrt(vecAB[0] ** 2 +
                                vecAB[1] ** 2 +
                                vecAB[2] ** 2)
        magnitudeCD = math.sqrt(vecCD[0] ** 2 +
                                vecCD[1] ** 2 +
                                vecCD[2] ** 2)

        # Compute cos(theta), where theta is the angle between the two vectors
        cosTheta = dotProduct / (magnitudeAB * magnitudeCD)

        # Return the angle
        return math.acos(cosTheta)

    def defineVector(self, coordsA, coordsB):
        '''
        Calculate a vector given two sets of coordinates
        Return the vector with its X, Y, Z coordinates
        '''

        vector = []

        # Calculating the vectors from the 3 points
        vector.append(coordsA[0] - coordsB[0])
        vector.append(coordsA[1] - coordsB[1])
        vector.append(coordsA[2] - coordsB[2])

        return vector

    def printFprint(self):
        '''
        Print the fprint, excluding the residue that don't have any interaction
        with the ligand
        '''

        if len(self.fprint) == 0:
            print("There is no fprint calculated. Generate fprints first.")
        else:
            print(",".join(self.fprint))

    def getFprint(self):
        '''
        Return the fprint, excluding residues that don't have any interaction
        with the ligand
        '''

        if len(self.fprint) == 0:
            print("There is no fprint calculated. Generate fprints first")
        else:
            return self.fprint

    def printFprintConsensus(self):
        '''
        Print the consensus fprint
        '''
        if len(self.fprintConsensus) == 0:
            print("There is no conensus fprint calculated, generate if first")
        else:
            print(",".join(self.fprintConsensus))

    def getFprintConsensus(self):
        '''
        Return the consensus fprint
        '''

        if len(self.fprintConsensus) == 0:
            print("There is no consensus fprint calculated, generate if first")
        else:
            return self.fprintConsensus

    def generateFprint(self):
        '''
        Loop over the self.residues and call all the interaction types defined
        in self.fprintDef
        '''
        # Hydrophobe parameters
        hydrophDist = 4.5

        # H-bond parameters
        # The cutoff defining H-Bond angle : ideal bond is 180deg (pi)
        # We allow +- 45deg (pi/4)
        # So the angle has to be less or equal to 225deg
        #   (180deg + 45deg, pi + pi/4 = 5/4pi)
        # And more or equal to 135deg
        #   (180deg - 45deg, pi - pi/4 = 3/4pi)
        hbondDist = 3.5
        hbondAngleUp = math.pi * (5.0 / 4)
        hbondAngleLow = math.pi * (3.0 / 4)

        # Weak h-bond parameters
        # The cutoff defining weak H-Bond angle :
        # Ideal bond is 90deg (pi),
        # And we allow +- 60deg (pi/3)
        # so the angle has to be <= to 150deg
        #   (90deg + 60deg, pi/2 + pi/3 = 5/6pi)
        # more or equal to 30deg
        #   (90deg - 60deg, pi/2 - pi/3 = pi/6)
        wkHbondDist = 2.8
        wkHbondAngleUp = math.pi * (5.0 / 6)
        wkHbondAngleLow = math.pi / 6.0

        # Cation Pi parameters
        catAniDist = 4.0

        # Aromatic parameters
        # The angle between rings is not checked, face-to-edge and face-to-face
        # interactions are grouped together
        aromDist = 5.0

        # Cation Pi parameters
        # Optimal angle is 180deg (pi) and
        # we allow a +- of pi/6 (30 deg)
        # The angle has to be between 7/6pi (210deg)
        # and 5/6pi (150 deg)
        catPiDist = 4.0
        catPiAngleUp = math.pi * (5.0 / 6)
        catPiAngleLow = math.pi / 6.0

        # Acceptor-Metal parameters
        accMetDist = 2.8

        # And setup: create fictive atoms that will represent the ring centers
        centerMol = oechem.OEGraphMol()
        resRingCenter = centerMol.NewAtom(0)
        ligRingCenter = centerMol.NewAtom(0)

        # Get a sorted list of the residues
        sortedKeys = self.molComplex.getResidues()

        # Check if there is a ligand: if there is no ligand, fill the
        # fingerprint with a fprint of zeros.
        if self.molComplex.getLigand() is None:
            for resTitle in sortedKeys:
                fprintSize = len(self.fprintDef.replace("0", ""))
                self.fprint.append(fprintSize * "0")
            return
        # If there is a ligand, proceed as usual, and calculate the fprint
        else:
            [ligand, ligRings] = self.molComplex.getLigand()

        for resTitle in sortedKeys:

            # Grab the residue from the dictionary
            residue, resRings, isConsensus = self.molComplex.residues[resTitle]

            # The fprint for this residue will be stored here, and appended
            # to the self.residues dictionary at the end of the loop
            resFprint = ""

            #-------------------------#
            # hydrophobe x hydrophobe #
            #-------------------------#
            if self.fprintDef[0] == "1":
                hydroph = self.check_inter(residue, ligand,
                                           self.ssHydrophobe, self.ssHydrophobe,
                                           hydrophDist)
                resFprint += hydroph

            #------------------------------#
            # donor (res) x acceptor (lig) #
            #------------------------------#
            if self.fprintDef[1] == "1":
                donAcc = self.check_hbond(residue, ligand,
                                          self.ssDonor, self.ssAcceptor,
                                          hbondDist, hbondAngleUp, hbondAngleLow)
                resFprint += donAcc

            #------------------------------#
            # donor (lig) x acceptor (res) #
            #------------------------------#
            if self.fprintDef[2] == "1":
                accDon = self.check_hbond(ligand, residue,
                                          self.ssDonor, self.ssAcceptor,
                                          hbondDist, hbondAngleUp, hbondAngleLow)
                resFprint += accDon

            #--------------------------------#
            # wkDon (res) x acc (lig)        #
            #   or wkDon (res) x wkAcc (lig) #
            #   or don (res) x wkAcc (lig)   #
            #--------------------------------#

            # Set wkDonAcc to 0, and check if each of the 3 interaction types
            # are found. Don't check the next interactions if one interaction
            # was already found

            if self.fprintDef[3] == "1":
                wkDonAcc = "0"

                # wkDon (res) x acc (lig)
                wkDonAcc = self.check_hbond(residue, ligand,
                                            self.ssWkDonor, self.ssAcceptor,
                                            wkHbondDist,
                                            wkHbondAngleUp, wkHbondAngleLow)
                if wkDonAcc == "0":
                    # wkDon (res) x wkAcc (lig)
                    wkDonAcc = self.check_hbond(residue, ligand,
                                                self.ssWkDonor, self.ssWkAcceptor,
                                                wkHbondDist,
                                                wkHbondAngleUp, wkHbondAngleLow)
                if wkDonAcc == "0":
                    # don (res) x wkAcc (lig)
                    wkDonAcc = self.check_hbond(residue, ligand,
                                                self.ssDonor, self.ssWkAcceptor,
                                                wkHbondDist,
                                                wkHbondAngleUp, wkHbondAngleLow)
                resFprint += wkDonAcc

            #--------------------------------#
            #  don (lig) x wkAcc (res)       #
            #   or wkDon (lig) x wkAcc (res) #
            #   or wkDon (lig) x acc (res)   #
            #--------------------------------#

            # Set wkDonAcc to 0, and check if each of the 3 interaction types
            # are found. Don't check the next interactions if one interaction
            # was already found

            if self.fprintDef[4] == "1":
                wkDonAcc = "0"

                # wkDon (lig) x acc (res)
                wkDonAcc = self.check_hbond(ligand, residue,
                                            self.ssDonor, self.ssWkAcceptor,
                                            wkHbondDist,
                                            wkHbondAngleUp, wkHbondAngleLow)
                if wkDonAcc == "0":
                    # wkDon (lig) x wkAcc (res)
                    wkDonAcc = self.check_hbond(ligand, residue,
                                                self.ssWkDonor, self.ssWkAcceptor,
                                                wkHbondDist,
                                                wkHbondAngleUp, wkHbondAngleLow)
                if wkDonAcc == "0":
                    # don (lig) x wkAcc (res)
                    wkDonAcc = self.check_hbond(ligand, residue,
                                                self.ssWkDonor, self.ssAcceptor,
                                                wkHbondDist,
                                                wkHbondAngleUp, wkHbondAngleLow)
                resFprint += wkDonAcc

            #----------------------------#
            # Cation (res) x Anion (lig) #
            #----------------------------#
            if self.fprintDef[5] == "1":
                catAni = self.check_inter(residue, ligand,
                                          self.ssCation, self.ssAnion,
                                          catAniDist)
                resFprint += catAni

            #----------------------------#
            # Anion (res) x Cation (lig) #
            #----------------------------#
            if self.fprintDef[6] == "1":
                aniCat = self.check_inter(residue, ligand,
                                          self.ssAnion, self.ssCation,
                                          catAniDist)
                resFprint += aniCat

            #----------------------------------------------------------#
            # Aromatic face2face and face2edge res x lig AND lig x res #
            #----------------------------------------------------------#
            if self.fprintDef[7] == "1":
                arom = self.check_arom(resRings, ligRings, centerMol,
                                       resRingCenter, ligRingCenter, aromDist)
                resFprint += arom

            #-------------------------#
            # Cation (res) x Pi (lig) #
            #-------------------------#
            if self.fprintDef[8] == "1":
                catPi = self.check_catPi(residue, self.ssCation,
                                         ligRings, centerMol, ligRingCenter,
                                         catPiDist, catPiAngleUp, catPiAngleLow)
                resFprint += catPi

            #-------------------------#
            # Pi (res) x Cation (lig) #
            #-------------------------#
            if self.fprintDef[9] == "1":
                piCat = self.check_catPi(ligand, self.ssCation,
                                         resRings, centerMol, resRingCenter,
                                         catPiDist, catPiAngleUp, catPiAngleLow)
                resFprint += piCat

            #------------------------------#
            # Acceptor (res) x Metal (lig) #
            #------------------------------#
            if self.fprintDef[10] == "1":
                accMet = self.check_inter(residue, ligand,
                                          self.ssAcceptor, self.ssMetal,
                                          accMetDist)
                resFprint += accMet

            self.fprint.append(resFprint)

    def check_inter(self, molA, molB,
                    patternA, patternB,
                    distCutoff):
        '''
        Function looping over the atoms of a residue and the ligand
        Depending on the pattern arguments, it will check for a type of
        interaction, and return 1 as soon as it found it, otherwise return 0
        '''

        # Look for pattern in residue
        for matchA in patternA.Match(molA):
            for atomA in matchA.GetTargetAtoms():

                # Look for pattern in ligand
                for matchB in patternB.Match(molB):
                    for atomB in matchB.GetTargetAtoms():

                        # Check distance between those atoms
                        dist = oechem.OEGetDistance(molA,
                                                    atomA,
                                                    molB,
                                                    atomB)

                        if dist <= distCutoff:
                            return "1"

        # Interaction not found between this residue and the ligand
        return "0"

    def check_hbond(self, molA, molB,
                    patternA, patternB,
                    distCutoff, angleCutUp, angleCutLow):
        '''
        Function checking hbond interaction between a residue and the ligand
        Checks for angle between:
            donor electronegative atom - donor proton - acceptor
        '''

        # Look for pattern in residue
        for matchA in patternA.Match(molA):
            for atomA in matchA.GetTargetAtoms():
                if str(atomA).split()[1] == 'H':
                    protonA = atomA

                    # Look for pattern in ligand
                    for matchB in patternB.Match(molB):
                        for acceptorB in matchB.GetTargetAtoms():

                            # Check distance between those atoms
                            dist = oechem.OEGetDistance(molA,
                                                        protonA,
                                                        molB,
                                                        acceptorB)

                            # if the distance meets the cutoff, check angle
                            if dist <= distCutoff:
                                #print dist, protonA, acceptorB, \
                                #      molA.GetTitle(), molB.GetTitle()
                                # Get the electronegative atom linked to
                                # the proton
                                for donorA in protonA.GetAtoms():
                                    pass
                                #print donorA
                                # Check angle: donorA, protonA, acceptorB
                                angle = oechem.OEGetAngle(molA, donorA,
                                                          molA, protonA,
                                                          molB, acceptorB)
                                #print angle, angleCutUp, angleCutLow
                                # Compare angle to cutoffs
                                if angle <= angleCutUp and \
                                        angle >= angleCutLow:
                                    return "1"

        # Interaction not found between this residue and the ligand
        return "0"

    def check_arom(self, resRings, ligRings, centerMol,
                   resRingCenter, ligRingCenter, aromDist):
        '''
        Checks for an interaction between the set of rings of a residue and a
        ligand. Uses temporary center atoms to make the distance comparisons
        '''

        # Loop over the residue rings
        for resRing in resRings:
            # For each ring, define its center using a temp atom
            centerMol.SetCoords(resRingCenter, resRing[0])

            # Loop over the ligand's rings
            for ligRing in ligRings:
                # Get the center of the ligand ring
                centerMol.SetCoords(ligRingCenter, ligRing[0])

                # Get the distance
                dist = oechem.OEGetDistance(centerMol, resRingCenter,
                                            centerMol, ligRingCenter)

                # If under cutoff, return 1
                if dist <= aromDist:
                    return "1"

        # If no aromatic interaction was found, return 0
        return "0"

    def check_catPi(self, mol, catPattern,
                    rings, centerMol, ringCenter,
                    catPiDist, catPiAngleUp, catPiAngleLow):
        '''
        Check for the cation-Pi interactions, loops over the rings of one
        molecule, and loops over the cations of the other
        '''

        # Loop over the rings
        for ring in rings:
            # Get the center of that ring
            centerMol.SetCoords(ringCenter, ring[0])

            # Loop over the cations of this 'molecule'
            for catMatch in catPattern.Match(mol):
                for cation in catMatch.GetTargetAtoms():

                    # Check distance between center of ring and cation
                    dist = oechem.OEGetDistance(centerMol, ringCenter,
                                                mol, cation)

                    #print cation, mol.GetTitle(), dist
                    # Check for distance threshold
                    if dist <= catPiDist:

                        # Calculate the vector between the ring's
                        # center and the cation (residue)
                        ringCenterCoords = centerMol.GetCoords(ringCenter)
                        cationCoords = mol.GetCoords(cation)
                        vecRingCation = self.defineVector(ringCenterCoords,
                                                          cationCoords)

                        # check the angle between the ring's normal,
                        # and the vector between the ring's center
                        # and the cation
                        angle = self.vectorAngle(vecRingCation, ring[2])

                        if angle <= catPiAngleUp or angle >= catPiAngleLow:
                            return "1"

        # If the interaction was not found, return 0
        return "0"
