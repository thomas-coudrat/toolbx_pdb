#!/usr/bin/env python

import confEnsemble

ens = confEnsemble.ConfEnsemble(
    "/home/thomas/Documents/FPRINT/fprintDataset")

ens.addTemplate("/home/thomas/Documents/FPRINT/template_08.pdb")

ens.makeComplexes()
ens.makeFprints()
ens.makeConsensusSeq()
ens.makeTanimoto("template_08.pdb")

ens.printDendrogram('jaccard')
#ens.printDendrogram('rogerstanimoto')

#ens.printFprints()
ens.printFprintsConsensus()

ens.makePCA("tanimoto")
ens.plotPCA("tanimoto", dim=2)

ens.plotPCA("tanimoto", dim=3)
