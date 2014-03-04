#!/usr/bin/env python

#------------------------------------------------------------------------------
#
# Work on a conformation ensemble in .pdb format and superimpose them all onto
# a user chosen .pdb file
#
# Thomas Coudrat, February 2014
#
#------------------------------------------------------------------------------

import os
import sys
import glob
import shutil
import argparse
from os.path import basename
from subprocess import check_output, STDOUT, CalledProcessError


def main():
    """
    Run the superimpose script
    """

    # Get the paths corresponding to the platform where this is executed
    icm, script = getPaths()
    # Get arguments
    pdbDir, templatePath = parsing()

    pdbPaths = glob.glob(pdbDir)

    print "Superimposing " + str(len(pdbPaths)) + " confs from dir: " + pdbDir
    print "Template used for superimposition: " + templatePath

    superimpose(templatePath, pdbPaths, icm, script)


def getPaths():
    """
    Figure the paths to the ICM executable and the ICM superimpose script
    """
    scriptDesktop = "/home/thomas/Copy/Tools/pdb_scripts/super.icm"
    icmDesktop = "/usr/icm-3.7-3b/icm64"
    scriptBarcoo = "/vlsci/VR0024/tcoudrat/Scripts/pdb_scripts/super.icm"
    icmBarcoo = "/vlsci/VR0024/tcoudrat/bin/icm-3.7-3b/icm64"
    #icm = "~/Install/ICM-3.7.3b/icm64"

    if os.path.exists(scriptDesktop):
        script = scriptDesktop
        icm = icmDesktop
    elif os.path.exists(scriptBarcoo):
        script = scriptBarcoo
        icm = icmBarcoo
    else:
        print "Error with ICM executable or ICM script"
        sys.exit()

    return icm, script


def parsing():
    """
    Parse arguments
    """
    descr = "Superimpose all .pdb conformations onto the provided template"
    descr_pdbDir = "Directory containing conformation ensemble"
    descr_templatePath = "Template onto which superimpose the conformation " \
        "ensemble"
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument("pdbDir", help=descr_pdbDir)
    parser.add_argument("templatePath", help=descr_templatePath)
    args = parser.parse_args()
    pdbDir = args.pdbDir
    templatePath = args.templatePath

    return pdbDir, templatePath


def superimpose(templatePath, pdbPaths, icm, script):
    """
    Get the template .pdb path as well as the paths to all the .pdb
    conformations
    Cycle through and for each of them run the ICM superimpose script.
    Make a temp copy of this script, and use sed to change for each execution
    """

    # Get the ICM name to be of the template
    templateName = basename(templatePath).replace(".pdb", "").replace("-", "_")

    for pdbPath in pdbPaths:
        pdbName = basename(pdbPath).replace(".pdb", "")
        # This is to reflect the fact that ICM turns hyphens to underscores
        # when a file is loaded in
        pdbName = pdbName.replace("-", "_")

        # Copy the superimpose ICM script to a temp location
        shutil.copy(script, "./temp.icm")

        # Make the sed changes
        os.system("sed -e 's|PDB_PATH_1|" + templatePath + "|g' ./temp.icm -i")
        os.system("sed -e 's|PDB_PATH_2|" + pdbPath + "|g' ./temp.icm -i")
        os.system("sed -e 's|PDB_NAME_1|" + templateName + "|g' ./temp.icm -i")
        os.system("sed -e 's|PDB_NAME_2|" + pdbName + "|g' ./temp.icm -i")

        print pdbName
        # Execute the temp script
        try:
            check_output(icm + " -s ./temp.icm", stderr=STDOUT, shell=True)
        except CalledProcessError, e:
            print "\n Error during the ICM script execution"
            print e.output
            sys.exit()

        # Remove the temp script
        os.remove("./temp.icm")


if __name__ == "__main__":
    main()
