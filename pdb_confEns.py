#!/usr/bin/env python

import confEnsemble
import argparse
import os


def main():
    """
    Run this script
    """

    ensDir, templatePath, top, dendro, pca, pcaLabel = parsing()

    ens = confEnsemble.ConfEnsemble(ensDir, top)
    ens.addTemplate(templatePath)

    ens.makeComplexes()
    ens.makeFprints()
    ens.makeConsensusSeq()
    ens.makeTanimoto(os.path.basename(templatePath))

    if dendro:
        ens.printDendrogram('jaccard')
        #ens.printDendrogram('rogerstanimoto')

    #ens.printFprints()
    ens.printFprintsConsensus()

    if pca:
        ens.makePCA("tanimoto")
        ens.plotPCA("tanimoto", dim=2)

    if pcaLabel:
        ens.plotPCA("tanimoto", dim=2, labelType="all")

    #ens.plotPCA("tanimoto", dim=3)


def parsing():
    """
    Parse arguments and define help file
    """

    descr = "This script executes commands on a protein conformation ensemble"
    descr_ensDir = "Path to the directory containing the .pdb files"
    descr_templatePath = "Path to the template that will be used to do " \
        "tanimoto comparisons on the interaction fingerprints"
    descr_top = "Select only the top X conformations from that " \
        "directory (optional)"
    descr_dendro = "Print-out a dendrogram of the conformations IFPs"
    descr_pca = "Print-out a PCA graph of the binding pocket conformations"
    descr_pcaLabel = "Print-out a PCA graph of the binding pocket " \
        "conformations (including labels on all conformations)"
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument("ensDir", help=descr_ensDir)
    parser.add_argument("templatePath", help=descr_templatePath)
    parser.add_argument("--top", help=descr_top)
    parser.add_argument("-dendro", action="store_true", help=descr_dendro)
    parser.add_argument("-pca", action="store_true", help=descr_pca)
    parser.add_argument("-pcaLabel", action="store_true", help=descr_pcaLabel)
    args = parser.parse_args()
    # Assign each variable parsed
    ensDir = args.ensDir
    templatePath = args.templatePath
    if args.top:
        top = int(args.top)
    else:
        top = None
    dendro = args.dendro
    pca = args.pca
    pcaLabel = args.pcaLabel

    return ensDir, templatePath, top, dendro, pca, pcaLabel


if __name__ == "__main__":
    main()
