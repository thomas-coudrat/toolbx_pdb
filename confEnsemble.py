#!/usr/bin/env python

# https://github.com/thomas-coudrat/toolbx_pdb
# Thomas Coudrat <thomas.coudrat@gmail.com>

import glob
import molecular_complex
import os
import sys
import fingerprint
import numpy
from scipy import spatial
from scipy import cluster
import matplotlib.pyplot as plt
import matplotlib as mpl
import PCA
import numpy as np
import math
import json


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

        # Get matplotlib to save SVG text as text, not paths
        mpl.rcParams['svg.fonttype'] = 'none'

        self.confDir = confDir
        confPaths = glob.glob(self.confDir + "/*.pdb")
        self.conformations = {}
        self.templates = []
        self.PCA = None
        self.consensusRes = None

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
                confName.replace(".pdb", "").replace("_", " ").upper()
                self.conformations[confName] = {'path': confPath}

    def initPCA(self):
        """
        Initiate PCA object
        """
        # Create PCA instance
        self.PCA = PCA.Principal_component_analysis(self.conformations)

    def generateProtCoords(self, consensusResidues=True):
        '''
        Extract protein coordinates data. If consensusResidues is True
        (default), then calculate PCA on the consensus residue list. Otherwise
        use the whole residue sequence.
        '''
        print("Calculating PCA on the conformation ensemble")
        self.PCA.makePCAcoords(consensusResidues)

    def calculate_and_plotPCA(self, projName, dim,
                              confLabels=None, metric=None):
        """
        Print or write plots for the PCA data loaded. Run this after
        self.makePCA has been ran.
        """
        if self.PCA is not None:
            # Create a labels list, with the name of the template only
            sortedConfNames = sorted(self.conformations.keys())
            # This variable will store labels to be used on the PCA score
            labels = []
            # Loop over conformation names
            for confName in sortedConfNames:
                # If that conformation is the template, add it to the labels
                if confName in self.templates:
                    labels.append(confName)
                # Otherwise, check if labels were provided as argument
                elif confLabels is not None:
                    # Display every labels if "all" is found
                    if confLabels[0] == "all":
                        labels.append(confName)
                    # Otherwise, display only the labels in that list
                    elif confName in confLabels:
                        labels.append(confName)
                    # Empty sting if no label was requested for this structure
                    else:
                        labels.append("")
                # If no label was requested, all structures get an empty string
                else:
                    labels.append("")

            # Call the plotPCAfig method
            self.PCA.plotPCAfig(projName, metric, labels, dim, self.templates)
        else:
            print("Create a PCA object with the function \
                  confEnsemble.makePCA() before plotting PCA data")
            sys.exit()

    def addConformation(self, pdbPath):
        '''
        Add a template .pdb file to the self.conformations path list. This
        template can be a crystal structure for example, adding it to the
        conformation lists will enable the same treatment on all structures,
        and keeping track of the template structure names will help make
        comparisons (e.g.: tanimoto comparison of fprint)
        '''

        if os.path.isfile(pdbPath) and pdbPath[-4:] == ".pdb":
            pdbName = os.path.basename(pdbPath)
            pdbName.replace(".pdb", "").replace("_", " ").upper()
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

    def generateQueryIFP(self, queryPath, fprintDef=None):
        """
        Generate the query IFP list and return it. Populate all residue IFP with
        'off' bits (e.g. 00000), and supplied residues with their corresponding
        IFP pattern.
        """

        # If fprint definition isn't supplied, default to full length IFPs
        if fprintDef is None:
            fprintDef = "111111111111"
        # Get the count of IFP booleans per residues
        resIFPcount = fprintDef.count("1")

        # Load json IFP query into a dictionary
        if os.path.isfile(queryPath) and queryPath[-5:] == ".json":
            queryName = os.path.basename(queryPath).replace(".json", "")
            with open(queryPath, "r") as jsonFile:
                queryDict = json.load(jsonFile)
        else:
            print("The following is not a .json file" + queryPath)
            sys.exit()

        # Get consensus residues from the first conformation (they all contain
        # the same consensus residues)
        first_conf_key = list(self.conformations.keys())[0]
        first_conf = self.conformations[first_conf_key]
        residues_consensus = first_conf['complex'].getResiduesConsensus()

        queryIFP = []
        # Loop over consensus residues
        for resName in residues_consensus:
            # If this consensus residue name is found in the query dictionary
            if resName in queryDict.keys():
                resIFP = queryDict[resName]
                # Check that it has the correct residue IFP count, if it does,
                # add it to the queryIFP list
                if len(resIFP) == resIFPcount:
                    queryIFP.append(resIFP)
                # Otherwise abort with error message
                else:
                    print(f"IFP residue pattern of '{resName}' doesn't " \
                          f"match the IFP length per residue of {resIFPcount}")
                    sys.exit()
            else:
                # Add an IFP residue pattern with booleans all turned to off '0'
                queryIFP.append(resIFPcount * "0")

        # Output progress
        print("\nQuery IFP data and IFP")
        print(queryDict)
        print(",".join(residues_consensus))
        print(",".join(queryIFP))

        return queryIFP

    def getTemplateIFP(self, templatePath):
        """
        Return the IFP from the molecule complex located in templatePath
        """
        # Format the template path into its name
        templateName = os.path.basename(templatePath)
        templateName.replace(".pdb", "").replace("_", " ").upper()

        # Get the template fprint
        if templateName in self.conformations:
            templateDict = self.conformations[templateName]
            templateFprintList = templateDict['fprint'].getFprintConsensus()
        else:
            print("\nTemplate not found: {}\n".format(templateName))
            sys.exit()

        return templateFprintList, templateName


    def computeDistances(self, ifpList, ifpName, metric):
        '''
        Calculate distance value between the chosen template
        (templatePos) and each of the conformations stored in
        self.conformations. Distance values are then stored in the
        corresponding dictionaries, under their corresponding distance name.
        '''

        # Go through all the conformations stored in self.conformations
        # Calculate the tanimoto coef between the chosen template and each
        # of their fprints, then store that coef in the conformationDict
        # with a new 'tanimoto' key
        for confName in sorted(self.conformations.keys()):
            # Get the dictionary of that conf
            conformationDict = self.conformations[confName]
            # Get the fprint list
            fprintList = conformationDict['fprint'].getFprintConsensus()

            # Compare two fprints using the tanimoto coefficient
            if metric == "tanimoto":
                tanimotoCoef = self.tanimoto(ifpList,
                                             fprintList)
                # Store this coef
                conformationDict['tanimoto'] = tanimotoCoef
            # Calculate Jaccard distance between two fprints
            elif metric == "jaccard":
                template_vect = np.array(list("".join(ifpList))).astype(bool)
                fprint_vect = np.array(list("".join(fprintList))).astype(bool)
                jaccardDist = spatial.distance.jaccard(template_vect,
                                                       fprint_vect)
                # print(template_vect, templateName)
                # print(fprint_vect, confName)
                # print(jaccardDist, "\n")
                conformationDict['jaccard'] = jaccardDist

            #print(ifpName, confName, jaccardDist)

        print(f"\nIFP {metric} distance:")
        print(f"- template: {ifpName}")

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

    def printDendrogram(self, projName, metric, dendroThresh, confLabels):
        '''
        Print a dendrogram based on the conensus fprints stored in
        self.conformations
        '''
        # If no dendrogram threshold was provided (None), default to a 0.5
        # threshold
        if dendroThresh is None:
            dendroThresh = 0.5

        print("\nGenerating dendrogram of conformations based on IFPs")
        print("Dendrogram threshold = {}\n".format(dendroThresh))

        # Temporarly store all the fprintCharStrings in this list
        fprintSpacedList = []
        confNamesList = []

        # Get all the fprints stored in self.conformations and modify
        # conformation names, ready for the figure
        for confName in sorted(self.conformations.keys()):
            # --------------------
            # Conformation names
            # --------------------
            # Prettify new name
            newConfName = confName.replace("_", " ").replace(".pdb", "").upper()
            # If conformation labels were provided, label them with *
            if confLabels is not None:
                if confName in confLabels:
                    newConfName = newConfName + " *"
            # Add all conformation "newNames" to the list
            confNamesList.append(newConfName)

            # -------------
            # Fingerprint
            # -------------
            # Get the dictionary of that conf
            conformationDict = self.conformations[confName]
            # Get the fprint list
            fprintList = conformationDict['fprint'].getFprintConsensus()
            # Create a fprint string
            fprintString = "".join(fprintList)
            # Insert a space between each of the bits in that new string
            # fprintSpacedList.append(list(fprintString).astype(bool))
            fprintCharList = list(fprintString)
            fprintSpaced = " ".join(fprintCharList)
            fprintSpacedList.append(fprintSpaced)

        # Create a concatenation of the spaced fprints
        fprintConcatenate = " ; ".join(fprintSpacedList)

        # ---------------------------
        # Creation of the dendrogram
        # ---------------------------

        # Create an array with the fprint concatenate
        mat = numpy.matrix(fprintConcatenate)
        X = numpy.asarray(mat)  # .astype(bool)

        # calculate the tanimoto coef into a matrix
        Y = spatial.distance.pdist(X, metric)

        # Building the hierachical clustering using the nearest point algo
        Z = cluster.hierarchy.linkage(Y, 'complete')
        # Clustering options are: single, complete, weighted, average

        # Setting parameters that cannot be modified via the SciPy API
        # This could be improved by using the matplotlib API via axes.
        mpl.rcParams['lines.linewidth'] = 5
        # Create the figure + axes
        fig = plt.figure()
        ax = fig.add_subplot(111)

        # Representation of the heirachical clustering in the form of a tree
        dendro = cluster.hierarchy.dendrogram(Z, ax=ax,
                                              color_threshold=dendroThresh,
                                              labels=confNamesList,
                                              leaf_font_size=15,
                                              orientation='left',
                                              above_threshold_color='black')

        # Add a vertical line for the threshold
        plt.axvline(dendroThresh, linewidth=3,
                    color='grey', linestyle="dashed")
        # plt.text(s=str(dendroThresh), x=dendroThresh, y=0.25,
        #         color="grey", fontsize=15, rotation="vertical")
        # plt.xticks(list(plt.xticks()[0]) + [dendroThresh])
        # plt.xticks([dendroThresh], [str(dendroThresh)], color="grey")

        # Changing figure style
        ax.tick_params(axis="x", which="major", labelsize=20)
        ax.spines['top'].set_color('none')
        ax.spines['left'].set_color('none')
        ax.spines['right'].set_color('none')
        ax.xaxis.set_ticks_position('bottom')
        if metric == "jaccard":
            ax.set_xlabel('Jaccard distance', fontsize=15)

        # Save the figure in svg format and png for quick visualization
        plt.savefig(projName + "_Dendro.svg", bbox_inches="tight")
        plt.savefig(projName + "_Dendro.png", bbox_inches="tight")

        # Resetting parameters
        mpl.rcParams['lines.linewidth'] = 1

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

    def makeFprints(self, fprintDef=None):
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
        Get the combination sequence of residues in all
        self.conformations that have a contact with the ligand
        '''

        consensusTemp = []

        # Extract a temporary list of all residues that make a contact with
        # the ligand, in any complex. Check by residueNumber_residueName, thus
        # storing only a single instance of each.
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
                res_num = res.split("_")[0]
                if res_num in [x.split("_")[0] for x in consensusTemp]:
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
        self.consensusRes = ",".join(consensusTemp)

    def printFprints(self):
        '''
        Print the fprints stored in self.conformations
        '''

        # Get the template fprints
        # for templateName in self.templates

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

    def plotFprints(self, projName, pdb_dir, fprintDef,
                    templatePath, additionalPaths):
        """
        Plot a graph representation of interaction fingerprints
        """
        # If no custom fingerprint was provided (None), then default to using
        # full featured interaction fingerprints
        if fprintDef is None:
            fprintDef = "111111111111"

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
        aa_codes = {"ALA": "A", "GLY": "G", "ILE": "I", "LEU": "L", "PRO": "P",
                    "VAL": "V", "PHE": "F", "TRP": "W", "TYR": "Y", "ASP": "D",
                    "GLU": "E", "ARG": "R", "HIS": "H", "LYS": "K", "SER": "S",
                    "THR": "T", "CYS": "C", "MET": "M", "ASN": "N", "GLN": "Q"}

        # Color palette for the current fprint definition
        colors_customFprint = []
        for bit, col in zip(fprintDef, colors_fullFprint):
            if bit == "1":
                colors_customFprint.append(col)

        # Get the conformations fprints
        fprintList = []
        # Initialise position lists and template positions
        addPathPos = []
        complPos = []
        templPathPos = False
        # Get the conformation list in descending order, when building the
        # figure that gets them ordered in ascending order
        sortedConfNames = sorted(self.conformations.keys(), reverse=True)
        # Loop over the conformations to get the data
        for i, confName in enumerate(sortedConfNames):
            # Get the dictionary of that conf
            conformationDict = self.conformations[confName]
            # Get the fprint list for "additional conformation(s)"
            confPath = conformationDict["path"]
            # Get the fprint list from other complexes
            fprintList.append(conformationDict['fprint'].getFprintConsensus())
            # Get the residue list
            resList = conformationDict['complex'].getResiduesConsensus()

            # Note the position of these for later renumbering
            if additionalPaths and confPath in additionalPaths:
                addPathPos.append(i)
            elif confPath == templatePath:
                templPathPos = i
            else:
                complPos.append(i)

        # Re-order the both the fprintList and the sortedConfNames to put the
        # additionalConformations first, followed by the template conformation
        # and finally
        if templPathPos:
            new_order_list = complPos + [templPathPos] + addPathPos
        else:
            new_order_list = complPos + addPathPos
        fprintList = [fprintList[i] for i in new_order_list]
        sortedConfNames = [sortedConfNames[i] for i in new_order_list]

        # Get the size of this custom fprint
        fp_length = len(fprintDef.replace("0", ""))
        # Get the number of residues contained in the fprint
        fp = fprintList[0]
        res_number = len(fp)
        # Get number of conformations
        confCount = len(self.conformations.keys())

        # Create the figure
        dpiVal = 800
        fig = plt.figure(figsize=(res_number/2, confCount), dpi=dpiVal)
        ax = fig.add_subplot(111)

        # Remove all but the squares
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.tick_params(left=False, top=False, bottom=False, right=False,
                       labelbottom=False, labelleft=False)

        # Create a set of spacers to be used to build the figure
        spacerX = 0.7
        # funky spacer calculation so that it gets calculated as a function of
        # the number of conformations
        spacerY = (1 - math.log(confCount, 10)) / 5000
        # Spacer multiplier in Y axis
        # spacerY = 0.00004
        # Generate default fprint X values. They get modified in the loop
        x_pos = []
        spacerCount = 0
        for i, bit in enumerate(np.arange(0, fp_length)):
            x_pos.append(spacerCount)
            spacerCount = spacerCount + spacerX
        x_pos = np.array(x_pos)

        for y, (fp, conf_name) in enumerate(zip(fprintList, sortedConfNames)):
            # Generate default fprint Y values positions
            y_pos = np.array([y*spacerY] * fp_length)

            # Write conformation names in front of each fprint scatter
            ax.text(-30, y_pos[0], conf_name.replace(".pdb", ""), size=4)

            # Loop over each fprint section representing a residue
            for i, (fp_segment, resName) in enumerate(zip(fp, resList)):

                # Create a set of colors corresponding to "on" bits of that
                # fp_segment and use colors defined above for the custom fprint
                colors = ["white"] * fp_length
                for j, bit in enumerate(fp_segment):
                    if bit == "1":
                        colors[j] = colors_customFprint[j]

                # Plot each fprint and add an offset of "fprint length" + 1
                # to include a gap between fp_segment
                ax.scatter(x_pos + (i * (fp_length + spacerX)),
                           y_pos,
                           s=5,
                           marker="s",
                           c=colors,
                           linewidths=0.1,
                           edgecolors="black")

                # Write residues only once (arbitrairly in first loop turn)
                if y == 0:
                    [resNumber, threeLetterCode] = resName.split("_")
                    # Get the single letter code for this residue
                    oneLetterCode = aa_codes[threeLetterCode]
                    # Change the number by removing the leading "0"s
                    if resNumber[0:2] == "00":
                        resNumber = resNumber[2]
                    elif resNumber[0:1] == "0":
                        resNumber = resNumber[1:3]
                    # Use this new residue name
                    new_resName = oneLetterCode + resNumber
                    ax.text(i * (fp_length + spacerX),
                            confCount * spacerY,
                            new_resName, size=5,  # rotation=90,
                            verticalalignment="bottom",
                            # horizontalalignment="right",
                            family="monospace")

        # Set limits on the figure
        ax.set_xlim([-20, len(fp) * (fp_length + spacerX)])
        # ax.set_ylim([-1 * spacerY, confCount * spacerY])

        # Save the figure in svg format and png for quick visualization
        plt.savefig(projName + "_IFP.svg", bbox_inches="tight")
        plt.savefig(projName + "_IFP.png", bbox_inches="tight", dpi=dpiVal)

    def printFprintsConsensus(self):
        '''
        Print the consensus fprints in self.conformations
        '''
        for confName in sorted(self.conformations):

            # Get the dictionary of that conf
            conformationDict = self.conformations[confName]
            # Get the residue list
            resList = conformationDict['complex'].getResiduesConsensus()
            # Get the ligand name in the complex
            ligObjs = conformationDict['complex'].getLigand()
            # Get the fprint list
            fprintList = conformationDict['fprint'].getFprintConsensus()

            # Print the data
            print(confName.replace(" ", "_"), ligObjs[0].GetTitle())
            print(",".join(resList))
            print(",".join(fprintList))

    def csvFprintsConsensus(self, projName, distance=False):
        '''
        Print the consensus fprints in self.conformations
        '''

        # Write CSV header, including consensus residue sequence
        with open("./" + projName + "_IFP.csv", "w") as f:
            if distance:
                f.write("Receptor,Ligand,IFPdist," + self.consensusRes + "\n")
            else:
                f.write("Receptor,Ligand," + self.consensusRes + "\n")

            for confName in sorted(self.conformations):
                # Get the dictionary of that conf
                conformationDict = self.conformations[confName]
                # Get the ligand name in the complex
                ligObjs = conformationDict['complex'].getLigand()
                # Get the fprint list
                fprintList = conformationDict['fprint'].getFprintConsensus()
                if distance:
                    # Get the jaccard distance between with the template IFP
                    IFPdist = conformationDict['jaccard']
                    # Write the data
                    f.write(confName.replace(" ", "_") + "," + ligObjs[0].GetTitle() + "," + str(IFPdist) + "," + ",".join(fprintList) + "\n")
                else:
                    # Write the data
                    f.write(confName.replace(" ", "_") + "," + ligObjs[0].GetTitle() + "," + ",".join(fprintList) + "\n")
