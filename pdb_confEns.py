#!/usr/bin/env python

# https://github.com/thomas-coudrat/toolbx_pdb
# Thomas Coudrat <thomas.coudrat@gmail.com>

import confEnsemble
import argparse
import os
import sys


def main():
    """
    Run this script
    """

    projName, ensDir, templatePath, additionalPaths, \
     top, dendro, dendroThresh, pca, pca3D, confLabels, \
     customFprint, ifp = parsing()

    # Write this command to a text file
    writeCommand(projName)

    # Create the conformation ensemble instance
    ens = confEnsemble.ConfEnsemble(ensDir, top)

    #-----------------------------------------
    # Generate interaction fingerprints (IFP)
    #-----------------------------------------

    # Add the template, this will add it to the list of complexes for the IFPs
    if templatePath:
        ens.addConformation(templatePath)
    if additionalPaths:
        for additionalPath in additionalPaths:
            ens.addConformation(additionalPath)
    # Generate molecular complexes (receptor/ligand)
    ens.makeComplexes()
    # Generate the interaction fingerprints
    ens.makeFprints(customFprint)
    # Generate a consensus sequence of residues among all conformations
    ens.makeConsensusSeq()
    # Save figure of interaction fingerprints representation
    if ifp:
        ens.plotFprints(projName, ensDir, customFprint,
                        templatePath, additionalPaths)
    # Print out IFPs
    # ens.printFprints()
    ens.printFprintsConsensus()
    ens.csvFprintsConsensus(projName)

    #------------------------------
    # Optional: Dendrogram and PCA
    #------------------------------

    # Display dendrogram
    if dendro:
        ens.printDendrogram(projName, 'jaccard', dendroThresh, confLabels)

    # Display PCA score (with or without labels)
    if pca:
        # Creates a PCA instance
        ens.initPCA()
        # Extract the protein coordinates on which PCA is calculated
        ens.generateProtCoords(consensusResidues=True)

        # Deal with conformation template for comparisons
        if templatePath:
            # Optionally compute distances between a template and conformations
            ens.computeDistances(templatePath, metric="jaccard")
            # Plot the PCA score
            ens.calculate_and_plotPCA(projName, dim=2,
                                      confLabels=confLabels,
                                      metric="jaccard")
            if pca3D:
                ens.calculate_and_plotPCA(projName, dim=3,
                                          confLabels=confLabels,
                                          metric="jaccard")
        else:
            ens.calculate_and_plotPCA(projName, dim=2,
                                      confLabels=confLabels)
            if pca3D:
                ens.calculate_and_plotPCA(projName, dim=3,
                                          confLabels=confLabels)
            #print("\nPCA not calculated: provide a template for " \
            #      "distance comparisons")
            #sys.exit()


def parsing():
    """
    Parse arguments and define help file
    """

    descr = "This script executes commands on a protein conformation ensemble"
    descr_projName = "Provide a project name. Format: string"
    descr_ensDir = "Path to the directory containing the .pdb files"
    descr_templatePath = "Path to the template that will be used to do " \
        "tanimoto comparisons on the interaction fingerprints"
    descr_additionalPaths = "Paths of additional conformations to be added " \
        "to the conformation ensemble analysis. " \
        "Format: path/to/conf1.pdb,path/to/conf2.pdb"
    descr_top = "Select only the top X conformations from that " \
        "directory (optional)"
    descr_dendro = "Print-out a dendrogram of the conformations IFPs"
    descr_dendroThresh = "Threshold to color the dendrogram. Value 0 < x < 1."
    descr_pca = "Print-out a PCA graph of the binding pocket conformations"
    descr_pca3D = "Print-out a PCA graph of the binding pocket conformations in 3D"
    descr_ifp = "Print-out an IFP diagram"
    descr_confLabel = "List pdb conformations to be identified in PCA and " \
        "Dendrogram plots. Format: 'conformation 1,conformation 4'"
    descr_customFprint = "Provide a custom interaction fingerprint " \
        "description of 11 bits (value 0 or 1), to inactivate or activate " \
        "of the following IFP descriptiors (in that order): " \
        "[hydrophobe x hydrophobe] " \
        "[donor (res) x acceptor (lig)] " \
        "[donor (lig) x acceptor (res)] " \
        "[wkDon (res) x acc (lig) or wkDon (res) x wkAcc (lig) or don (res) x wkAcc (lig)] " \
        "[don (lig) x wkAcc (res) or wkDon (lig) x wkAcc (res) or wkDon (lig) x acc (res)] " \
        "[Cation (res) x Anion (lig)] " \
        "[Anion (res) x Cation (lig)] " \
        "[Aromatic face2face and face2edge res x lig AND lig x res] " \
        "[Cation (res) x Pi (lig)] " \
        "[Pi (res) x Cation (lig)] " \
        "[Acceptor (res) x Metal (lig)]"
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument("projName", help=descr_projName)
    parser.add_argument("ensDir", help=descr_ensDir)
    parser.add_argument("-templatePath", help=descr_templatePath)
    parser.add_argument("-additionalPaths", help=descr_additionalPaths)
    parser.add_argument("--top", help=descr_top)
    parser.add_argument("-dendro", action="store_true", help=descr_dendro)
    parser.add_argument("--dendroThresh", help=descr_dendroThresh)
    parser.add_argument("-pca", action="store_true", help=descr_pca)
    parser.add_argument("-pca3D", action="store_true", help=descr_pca3D)
    parser.add_argument("-ifp", action="store_true", help=descr_ifp)
    parser.add_argument("--confLabels", help=descr_confLabel)
    parser.add_argument("-customFprint", help=descr_customFprint)
    args = parser.parse_args()

    #-----------------------------
    # Assign each variable parsed
    #-----------------------------

    projName = args.projName.replace(" ", "_")

    ensDir = args.ensDir

    templatePath = args.templatePath

    if args.additionalPaths:
        additionalPaths = args.additionalPaths.split(",")
        if not all([x[-4:] == (".pdb") for x in additionalPaths]):
            print("Additional conformation paths must be .pdb files. Exiting.")
            sys.exit()
    else:
        additionalPaths = None

    if args.top:
        top = int(args.top)
    else:
        top = None

    dendro = args.dendro

    if args.dendroThresh:
        # Make sure that the threshold value is a float
        try:
            dendroThresh = float(args.dendroThresh)
        except ValueError:
            print("Dendrogram threshold should be a float/int")
        # Make sure that the threshold value submitted is between 0 and 1
        if 0.0 <= dendroThresh <= 1.0:
            pass
        else:
            print("Dendrogram threshold must have numerical value 0 < x < 1")
            sys.exit()
    else:
        dendroThresh = None

    pca = args.pca

    pca3D = args.pca3D

    ifp = args.ifp

    if args.confLabels:
        confLabels = args.confLabels.split(",")
    else:
        confLabels = None

    if args.customFprint:
        customFprint = args.customFprint
        if len(customFprint) != 11:
            print("Input error: custom fprint definition has to be 11 bits long")
            sys.exit()
        if len(customFprint.replace("0", "").replace("1", "")) != 0:
            print("Input error: custom fprint definition uses only '0' and/or '1'")
            sys.exit()
    else:
        customFprint = None

    return projName, ensDir, templatePath, additionalPaths, \
        top, dendro, dendroThresh, pca, pca3D, confLabels, customFprint, ifp


def writeCommand(projName):
    """
    Write down the command that was used to exectute this script in a .sh
    file, at the location where the script is executed. Also write the
    current working directory at the time of execution
    """

    cwd = os.getcwd()
    logFile = open(projName + "_CMD.sh", "w")
    # Write the directory location: this is not executed upong sh call of
    # the thisFile.sh, but serves as information
    logFile.write(cwd + "\n")
    logFile.write(sys.argv[0].split("/")[-1] + " ")
    for arg in sys.argv[1:]:
        if len(arg) > 0:
            # Deal with argument options (starting with '-')
            if arg[0] == "-":
                logFile.write(arg + " ")
            # Do not add "'" on argument if it already has them
            elif arg[0] == "'" and arg[-1] == "'":
                logFile.write(arg + " ")
            # Add the "'" around each other argument
            else:
                logFile.write("'" + arg + "' ")
    logFile.close()


if __name__ == "__main__":
    main()
