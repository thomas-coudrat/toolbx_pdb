#!/usr/bin/env python

import confEnsemble
import argparse
import os


def main():
    """
    Run this script
    """

    ensDir, topX, templatePath = parsing()

    ens = confEnsemble.ConfEnsemble(ensDir, topX)
    ens.addTemplate(templatePath)

    ens.makeComplexes()
    ens.makeFprints()
    ens.makeConsensusSeq()
    ens.makeTanimoto(os.path.basename(templatePath))

    ens.printDendrogram('jaccard')
    #ens.printDendrogram('rogerstanimoto')

    #ens.printFprints()
    ens.printFprintsConsensus()

    ens.makePCA("tanimoto")
    ens.plotPCA("tanimoto", dim=2)
    ens.plotPCA("tanimoto", dim=2, labelType="all")
    #ens.plotPCA("tanimoto", dim=3)


def parsing():
    """
    Parse arguments and define help file
    """

    descr = "This script executes commands on a protein conformation ensemble"
    descr_ensDir = "Path to the directory containing the .pdb files"
    descr_topX = "Select only the top X conformations from that" \
        "directory (optional)"
    descr_templatePath = "Path to the template that will be used to do" \
        "tanimoto comparisons on the interaction fingerprints"
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument("ensDir", help=descr_ensDir)
    parser.add_argument("--topX", help=descr_topX)
    parser.add_argument("templatePath", help=descr_templatePath)
    args = parser.parse_args()
    ensDir = args.ensDir
    topX = int(args.topX)
    templatePath = args.templatePath

    return ensDir, topX, templatePath


if __name__ == "__main__":
    main()
