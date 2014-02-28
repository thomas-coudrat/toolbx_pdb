#!/usr/bin/env python

import molecular_complex
import fingerprint
import confEnsemble


pdbPath = "pdbs/ldm_03.pdb"

macro = molecular_complex.Complex(pdbPath)
#print macro
macro.makeComplex()
#macro.printResidues()
#macro.printLigand()

# Depreciated
#resList, fprintList = macro.getFprint()
#for r, f in zip(resList, fprintList):
#    print r, f

fprint = fingerprint.Fingerprint(macro)
#print fprint.fprint
fprint.generateFprint()
#fprint.printFprint("column")
#fprint.printFprint("string")
#resList, fprintList = fprint.getFprint()

#for r in resList:
#    print r

#macro.printResidues()

resList, fprintList = macro.getFprint()
for r, f in zip(resList, fprintList):
    print r, f
