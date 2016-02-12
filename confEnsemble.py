#!/usr/bin/env python

import glob
import molecular_complex
import os
import sys
import fingerprint
import numpy
from scipy import spatial
from scipy import cluster
import matplotlib.pyplot as plt
import PCA
import numpy as np


class ConfEnsemble:
    """
    Class holds a dictionary self.conformations that gets created with
    'conformation names' as a key.
    Each of the objects stored with these keys are dictionaries, that get
    initialized with the 'path' of the conformation.
    This conformation dictionary can hold those keys if created:
        - 'path' = string
        - 'complex' = molecular_complex.Complex object
        - 'fprint' = fingerprint.Fingerprint object
        - 'tanimoto' = float
    """

    def __init__(self, confDir, topX=None):
        '''
        Create a Conformations instance, with a directory path
        '''
        self.confDir = confDir
        confPaths = glob.glob(self.confDir + "/*.pdb")
        self.conformations = {}
        self.templates = []
        self.PCA = None

        # Check if any .pdb file were found
        if len(confPaths) == 0:
            print("This directory does not contains any .pdb files!")
        else:
            # Initialize the self.conformations dictionary, each key points to
            # a dictionary to store the various data related to each
            # conformation
            if topX:
                # Only the top X conformations, sorted by alphanumeric order,
                # will be considered if a topX value is provided
                confPaths = sorted(confPaths)
                confPaths = confPaths[0:topX]

            for confPath in confPaths:
                confName = os.path.basename(confPath)
                self.conformations[confName] = {'path': confPath}

    def makePCA(self, var_to_plot=None, residues=None):
        '''
        Initiate the creation of a PCA object, that will get the PCA for the
        conformation ensemble
        'residues' is a string containing the residues to be included in the
        PCA. Consider the whole sequence if residues=None
        '''
        if residues is not None:
            l_resList = residues.split(",")
            for i, res in enumerate(l_resList):
                l_resList[i] = int(res.split("_")[0])
        else:
            l_resList = None

        self.PCA = PCA.Principal_component_analysis(self.conformations,
                                                    l_resList)
        self.PCA.makePCAcoords()
        if var_to_plot is not None:
            self.PCA.makePCAvars(var_to_plot)

    def plotPCA(self, var_to_plot, dim, labelType="template"):
        """
        Print or write plots for the PCA data loaded. Run this after
        self.makePCA has been ran.
        """
        if self.PCA is not None:
            # Create a labels list, with the name of the template only
            sortedConfNames = sorted(self.conformations.keys())
            labels = []
            for confName in sortedConfNames:
                if labelType == "template":
                    if confName in self.templates:
                        labels.append(confName.replace(".pdb", ""))
                    else:
                        labels.append("")
                elif labelType == "all":
                    labels.append(confName.replace(".pdb", ""))
                else:
                    print("This is not a recognised 'labelType' option to " +
                          "plot the PCA coordinates:" + labelType)
                    sys.exit()
            # Call the plotPCAfig method
            self.PCA.plotPCAfig(var_to_plot, labels, dim)
        else:
            print("You first have to create a PCA object with the function")
            print("Use the function confEnsemble.makePCA()")

    def addTemplate(self, pdbPath):
        '''
        Add a template .pdb file to the self.conformations path list. This
        template can be a crystal structure for example, adding it to the
        conformation lists will enable the same treatment on all structures,
        and keeping track of the template structure names will help make
        comparisons (e.g.: tanimoto comparison of fprint)
        '''

        if os.path.isfile(pdbPath):
            pdbName = os.path.basename(pdbPath)
            self.conformations[pdbName] = {'path': pdbPath}
            self.templates.append(pdbName)
        else:
            print("The following is not a .pdb file" + pdbPath)

    def makeFakeValue(self):
        '''
        Add a fake value to each conformation dict in order to be able to
        use this class, when no coloring is displayed on the PCA graph
        '''

        for confName in self.conformations:
            conformationDict = self.conformations[confName]
            conformationDict['fake'] = 0

    def makeTanimoto(self, templateName):
        '''
        Calculate tanimoto coefficient between the chosen template
        (templatePos) and each of the conformations stored in
        self.conformations. The tanimoto coefficient is then stored in the
        corresponding dictionaries, with the key 'tanimoto'.
        '''

        # Get the template fprint
        if templateName in self.conformations:
            templateDict = self.conformations[templateName]
            templateFprintList = templateDict['fprint'].getFprintConsensus()

            # Go through all the conformations stored in self.conformations
            # Calculate the tanimoto coef between the chosen template and each
            # of their fprints, then store that coef in the conformationDict
            # with a new 'tanimoto' key
            for confName in self.conformations:

                # Get the dictionary of that conf
                conformationDict = self.conformations[confName]
                # Get the fprint list
                fprintList = conformationDict['fprint'].getFprintConsensus()

                # Calculate the Tanimoto coef between those two fprints
                tanimotoCoef = self.tanimoto(templateFprintList, fprintList)

                # Store this coef
                conformationDict['tanimoto'] = tanimotoCoef

            print("Tanimoto coefficients were calculated using template:" +
                  templateName)
        else:
            print("This template was not found for tanimoto calculations" +
                  templateName)

    def tanimoto(self, fprintListA, fprintListB):
        '''
        Return the tanimoto coefficient rounded to the 4th decimal, between
        the two fprints provided
        '''
        # initialize numerator and denominator of the tanimoto coefficient
        # between the profile and current fingerprint (fprint)
        num = 0
        denom = 0
        # Make the fprints strings
        fprintA = "".join(fprintListA)
        fprintB = "".join(fprintListB)

        # Loop over each bit in current fprint
        for k in range(len(fprintA)):
            num += int(fprintA[k]) * int(fprintB[k])
            denom += int(fprintA[k]) ** 2 + int(fprintB[k]) ** 2 - \
                int(fprintA[k]) * int(fprintB[k])

        # Calculate tanimoto only is denominator is not 0
        if denom == 0:
            tanimoto = 0
        else:
            tanimoto = num * 1.0 / denom

        return 1 - round(tanimoto, 4)

    def printDendrogram(self, metric):
        '''
        Print a dendrogram based on the conensus fprints stored in
        self.conformations
        '''

        # Temporarly store all the fprintCharStrings in this list
        fprintSpacedList = []
        confNamesList = []

        # Get all the fprints stored in self.conformations
        for confName in self.conformations:

            # Keep the conformation names in a list, to rename the dendrogram
            # later
            confNamesList.append(confName)

            # Get the dictionary of that conf
            conformationDict = self.conformations[confName]
            # Get the fprint list
            fprintList = conformationDict['fprint'].getFprintConsensus()
            # Create a fprint string
            fprintString = "".join(fprintList)

            # Insert a space between each of the bits in that new string
            fprintCharList = list(fprintString)
            fprintSpaced = " ".join(fprintCharList)
            fprintSpacedList.append(fprintSpaced)

        # Create a concatenation of the spaced fprints
        fprintConcatenate = " ; ".join(fprintSpacedList)

        # Creation of the dendrogram
        #---------------------------

        # Create an array with the fprint concatenate
        mat = numpy.matrix(fprintConcatenate)
        X = numpy.asarray(mat)

        # calculate the tanimoto coef into a matrix
        Y = spatial.distance.pdist(X, metric)
        # Here a custom tanimoto function could be used

        # Building the hierachical clustering using the nearest point algo
        Z = cluster.hierarchy.linkage(Y, 'complete')
        # other options are: single, complete, weighted, average

        # Representation of the heirachical clustering in the form of a tree
        plt.figure()
        dendro = cluster.hierarchy.dendrogram(Z,
                                              color_threshold=0.5,
                                              orientation='right')

        # Adding the names of the conformations to the dendrogram before
        # showing it
        for i, pos in enumerate(dendro['ivl']):
                dendro['ivl'][int(i)] = confNamesList[int(pos)]

        # Display the dendrogram
        plt.show()

    def plotSortedEnsemble(self, val):
        """
        Print the sorted ensemble of conformations in the ensemble, sorted
        along the X axis, with a value given as parameter 'val' to be displayed
        along the Y axis
        """

        fig = plt.figure()
        axes = fig.add_subplot(111)

        sortedKeys = self.conformations.keys()
        sortedKeys.sort()

        xData = []
        yData = []
        for i, confName in enumerate(sortedKeys):
            conformationDict = self.conformations[confName]
            value = conformationDict[val]
            print(confName + " " + value)
            xData.append(i)
            yData.append(value)

        axes.plot(xData, yData)

        plt.show()

    def makeComplexes(self):
        '''
        Initiate the creation of Complex objects, for all the conformations
        stored in self.conformations
        '''

        print("\nCreating molecular complex objects:")

        total = len(self.conformations)
        for i, confName in enumerate(self.conformations):

            # Get the dictionary associated with this conformation
            conformationDict = self.conformations[confName]
            # Get the path of the conformation from the dictionary
            confPath = conformationDict['path']
            # Create a molecular complex
            molComplex = molecular_complex.Complex(confPath)
            # Populate the molecular complex object
            molComplex.makeComplex()
            # store it in the dictionary
            conformationDict['complex'] = molComplex
            print("Mol-complex: " + str(i + 1) + "/" + str(total))

    def makeFprints(self, fprintDef=False):
        '''
        Initiate the creation of fprints for all the conformations stored in
        self.conformations dictionary
        '''

        print("\nGenerating ligand/protein interaction fingerprints:")

        total = len(self.conformations)
        for i, confName in enumerate(self.conformations):

            # Get the dictionary associated with this conformation
            conformationDict = self.conformations[confName]
            # Get the molecular complex of that conformation
            molComplex = conformationDict['complex']
            # Apply the fingerprint on that complex
            fprint = fingerprint.Fingerprint(molComplex, fprintDef)
            # Populate the molComplex object with the fprints
            fprint.generateFprint()
            # Store this fprint object into the corresponding confDict dict
            conformationDict['fprint'] = fprint
            print("fprint: " + str(i + 1) + "/" + str(total))

    def makeConsensusSeq(self):
        '''
        Get a consensus sequence of residues in all self.conformations that
        have a contact with the ligand
        '''

        consensusTemp = []

        for confName in self.conformations:
            # Get the dictionary of that conf
            conformationDict = self.conformations[confName]
            # Get the residue list
            resList = conformationDict['complex'].getResidues()
            # Get the fprint list
            fprintList = conformationDict['fprint'].getFprint()

            for res, fprint in zip(resList, fprintList):
                if "1" in fprint and res not in consensusTemp:
                    consensusTemp.append(res)

        # Not essential to have it sorted, but do it anyway
        consensusTemp.sort()

        # Get the length of a residue fprint, and create a residue blank fprint
        fprintLen = len(fprint)
        fprintBlank = fprintLen * "0"

        for confName in self.conformations:

            # get the conformation dict
            conformationDict = self.conformations[confName]
            # Get the molComplex object
            molComplex = conformationDict['complex']
            # Get the residue list
            resListSorted = molComplex.getResidues()
            # Get the fprint object
            fprintObject = conformationDict['fprint']
            # Get the fprint list
            fprintList = fprintObject.getFprint()

            fprintConsensus = []
            for res, fprint in zip(resListSorted, fprintList):
                if res in consensusTemp:
                    if "1" in fprint:
                        fprintConsensus.append(fprint)
                    else:
                        fprintConsensus.append(fprintBlank)
                else:
                    # Every time a residue is not part of the consensus seq
                    # add a False flag to the molComplex.residues dictionary
                    residue, resRings, isConsensus = molComplex.residues[res]
                    molComplex.residues[res] = [residue, resRings, False]

            # Store the new consensusFprint in the fprint object
            fprintObject.fprintConsensus = fprintConsensus

        consensusTemp.sort()
        consensusLen = str(len(consensusTemp))
        print("\nConsensus sequence of length " + consensusLen + " generated:")
        print(",".join(consensusTemp) + "\n")

    def printFprints(self):
        '''
        Print the fprints stored in self.conformations
        '''

        # Get the template fprints
        #for templateName in self.templates

        # Get the conformations fprints
        for confName in self.conformations:

            # Get the dictionary of that conf
            conformationDict = self.conformations[confName]
            # Get the residue list
            resList = conformationDict['complex'].getResidues()
            # Get the fprint list
            fprintList = conformationDict['fprint'].getFprint()

            # Print the data
            print(confName + " " + str(len(resList)))
            print(",".join(resList))
            print(",".join(fprintList) + "\n")

    def plotFprints(self, fprintDef="111111111111"):
        """
        Plot a graph representation of interaction fingerprints
        """

        # Color definition for the full fprint definition
        c_hydroph = "blue"
        c_hbond = "red"
        c_wkHbond = "orange"
        c_ionic = "purple"
        c_aromatic = "green"
        c_catPi = "pink"
        c_metal = "yellow"
        colors_fullFprint = [c_hydroph,
                             c_hbond,
                             c_hbond,
                             c_wkHbond,
                             c_wkHbond,
                             c_ionic,
                             c_ionic,
                             c_aromatic,
                             c_catPi,
                             c_catPi,
                             c_metal]

        # Color palette for the current fprint definition
        colors_customFprint = []
        for bit, col in zip(fprintDef, colors_fullFprint):
            if bit == "1":
                colors_customFprint.append(col)

        # Get the conformations fprints
        fprintList = []
        sortedConfNames = sorted(self.conformations.keys())
        for confName in sortedConfNames:
            # Get the dictionary of that conf
            conformationDict = self.conformations[confName]
            # Get the residue list
            resList = conformationDict['complex'].getResiduesConsensus()
            # Get the fprint list
            fprintList.append(conformationDict['fprint'].getFprintConsensus())

        # Get the size of this custom fprint
        fp_length = len(fprintDef.replace("0", ""))
        # Get the number of residues contained in the fprint
        fp = fprintList[0]
        res_number = len(fp)
        # Get number of conformations
        confCount = len(self.conformations.keys())

        # Create the figure
        dpiVal = 800
        fig = plt.figure(figsize=(res_number/2,confCount*1), dpi=dpiVal)
        ax = fig.add_subplot(111)

        # Remove all but the squares
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.tick_params(left="off", top="off", bottom="off", right="off",
                       labelbottom='off', labelleft='off')

        spacerX = 1
        spacerY = 0.0002
        for y, (fp, conf_name) in enumerate(zip(fprintList, sortedConfNames)):
            # Generate default fprint X and Y values positions
            # (X values get modified in the for loop)
            x_pos = np.arange(0, fp_length)
            y_pos = np.array([y*spacerY] * fp_length)

            # Write conformation names in front of each fprint scatter
            ax.text(-30, y_pos[y], conf_name.replace(".pdb",""), size=8)

            # Loop over each fprint section representing a residue
            for i, (fp_segment, resName) in enumerate(zip(fp, resList)):

                # Create a set of colors corresponding to "on" bits of that
                #fp_segment and use colors defined above for the custom fprint
                colors = ["white"] * fp_length
                for j, bit in enumerate(fp_segment):
                    if bit == "1":
                        colors[j] = colors_customFprint[j]

                # Plot each fprint and add an offset of "fprint length" + 1
                # to include a gap between fp_segment
                ax.scatter(x_pos + (i * (fp_length + spacerX)),
                           y_pos,
                           s=15,
                           marker="s",
                           c=colors,
                           linewidths=0.3)

                # Write residues only once (arbitrairly in first loop turn)
                if y == 0:
                    ax.text(i * (fp_length + spacerX) + (fp_length/2 - spacerX),
                            confCount * spacerY,
                            resName.replace("_", " "), size=8, rotation=90,
                            verticalalignment="bottom",
                            family="monospace")

        # Set limits on the figure
        ax.set_xlim([-20, len(fp) * (fp_length + spacerX)])
        #ax.set_ylim([-1 * spacerY, confCount * spacerY])

        # Save the figure in pdf format
        plt.savefig("ifp.pdf", bbox_inches="tight", format="pdf", dpi=dpiVal)

    def printFprintsConsensus(self):
        '''
        Print the consensus fprints in self.conformations
        '''
        for confName in self.conformations:

            # Get the dictionary of that conf
            conformationDict = self.conformations[confName]
            # Get the residue list
            resList = conformationDict['complex'].getResiduesConsensus()
            # Get the fprint list
            fprintList = conformationDict['fprint'].getFprintConsensus()

            # Print the data
            print(confName)
            #print(",".join(resList))
            print(",".join(fprintList))
