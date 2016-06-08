#!/usr/bin/env python

import confEnsemble
import argparse
import os
import sys


def main():
    """
    Run this script
    """

    projName, ensDir, templatePath, top, dendro, dendroThresh, pca, pcaLabels, customFprint = parsing()

    # Create the conformation ensemble instance
    ens = confEnsemble.ConfEnsemble(ensDir, top)

    #-----------------------------------------
    # Generate interaction fingerprints (IFP)
    #-----------------------------------------

    # Add the template, this will add it to the list of complexes for the IFPs
    if templatePath:
        ens.addTemplate(templatePath)
    # Generate molecular complexes (receptor/ligand)
    ens.makeComplexes()
    # Generate the interaction fingerprints
    ens.makeFprints(customFprint)
    # Generate a consensus sequence of residues among all conformations
    ens.makeConsensusSeq()
    # Save figure of interaction fingerprints representation
    if customFprint:
        ens.plotFprints(projName, ensDir, customFprint)
    else:
        ens.plotFprints(projName, ensDir)
    # Print out IFPs
    # ens.printFprints()
    ens.printFprintsConsensus()

    #------------------------------
    # Optional: Dendrogram and PCA
    #------------------------------

    # Display dendrogram
    if dendro:
        if dendroThresh:
            ens.printDendrogram(projName, 'jaccard', dendroThresh)
            # Other metric that could be used:
            #ens.printDendrogram('rogerstanimoto', dendroThresh)
        else:
            ens.printDendrogram(projName, 'jaccard')

    # Display PCA score (with or without labels)
    if pca:
        if templatePath:
            # Calculating tanimoto comparisons is required for the PCA score
            ens.computeDistances(os.path.basename(templatePath), "jaccard")
            ens.makePCA("jaccard")
            if pcaLabels:
                ens.plotPCA(projName, "jaccard", dim=2, pcaLabels=pcaLabels)
            else:
                ens.plotPCA(projName, "jaccard", dim=2)
        else:
            print("\nPCA not calculated: provide a template for " \
                  "tanimoto comparisons")
            sys.exit()

    writeCommand(projName)


def parsing():
    """
    Parse arguments and define help file
    """

    descr = "This script executes commands on a protein conformation ensemble"
    descr_projName = "Provide a project name. Format: string"
    descr_ensDir = "Path to the directory containing the .pdb files"
    descr_templatePath = "Path to the template that will be used to do " \
        "tanimoto comparisons on the interaction fingerprints"
    descr_top = "Select only the top X conformations from that " \
        "directory (optional)"
    descr_dendro = "Print-out a dendrogram of the conformations IFPs"
    descr_dendroThresh = "Threshold to color the dendrogram. Value 0 < x < 1."
    descr_pca = "Print-out a PCA graph of the binding pocket conformations"
    descr_pcaLabel = "List pdb conformations to be labelled in PCA plot. " \
        "Format: conformation1.pdb,conformation4.pdb"
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
    parser.add_argument("--top", help=descr_top)
    parser.add_argument("-dendro", action="store_true", help=descr_dendro)
    parser.add_argument("--dendroThresh", help=descr_dendroThresh)
    parser.add_argument("-pca", action="store_true", help=descr_pca)
    parser.add_argument("--pcaLabels", help=descr_pcaLabel)
    parser.add_argument("-customFprint", help=descr_customFprint)
    args = parser.parse_args()

    #-----------------------------
    # Assign each variable parsed
    #-----------------------------

    projName = args.projName.replace(" ", "_")

    ensDir = args.ensDir

    templatePath = args.templatePath

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
        dendroThresh = False

    pca = args.pca

    if args.pcaLabels:
        pcaLabels = args.pcaLabels.split(",")
        if not all([x[-4:] == (".pdb") for x in pcaLabels]) and not pcaLabels[0] == "all" :
            print("PCA label must be .pdb conformations. Exiting.")
            sys.exit()
    else:
        pcaLabels = None

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

    return projName, ensDir, templatePath, top, dendro, dendroThresh, pca, pcaLabels, customFprint


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
