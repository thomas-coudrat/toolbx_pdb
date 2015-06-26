#!/usr/bin/env python

# ------------------------------------------------------------------------------
#
# Work on a conformation ensemble in .pdb format and superimpose them all onto
# a user chosen .pdb file
#
# Thomas Coudrat, February 2014
# thomas.coudrat@gmail.com
#
# ------------------------------------------------------------------------------

import os
import sys
import glob
import shutil
import argparse
from os.path import basename
from subprocess import check_output, STDOUT, CalledProcessError
import socket
import json


def main():
    """
    Run the superimpose script
    """

    # Get arguments
    pdbDir, templatePath = parsing()

    # Get the paths corresponding to the platform where this is executed
    icm, script = getPaths()

    # Get paths of all .pdb files in this directory
    pdbPaths = sorted(glob.glob(pdbDir + "/*.pdb"))

    print
    print "Superimposing " + str(len(pdbPaths)) + " confs from dir: " + pdbDir
    print "Template used for superimposition: " + templatePath
    print

    superimpose(templatePath, pdbPaths, pdbDir, icm, script)


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


def getPaths():
    """
    Figure the paths to the ICM executable and the ICM superimpose script
    """

    # This Json file stores the ICM executable locations for each platform
    hostFiles_json = os.path.dirname(os.path.realpath(__file__)) + \
        "/hostFiles.json"

    # Read content of .json file
    with open(hostFiles_json, "r") as jsonFile:
        hostFiles_dict = json.load(jsonFile)

    # Get the hostname
    hostname = socket.gethostname()

    # Verify that the current hostname is present in the config file loaded from
    # json
    if hostname in hostFiles_dict.keys():
        icm = hostFiles_dict[hostname][0]
        script = hostFiles_dict[hostname][1]
    else:
        print "Error: hostname " + hostname + " filepaths not configured in " \
            " ./hostFiles.json"
        sys.exit()

    return icm, script


def superimpose(templatePath, pdbPaths, pdbDir, icm, script):
    """
    Get the template .pdb path as well as the paths to all the .pdb
    conformations
    Cycle through and for each of them run the ICM superimpose script.
    Make a temp copy of this script, and use sed to change for each execution
    """

    # Get the ICM name to be of the template
    templateName = basename(templatePath).replace(".pdb", "").replace("-", "_")
    # Temp directory name
    tempScript = "./" + os.path.dirname(pdbDir).lower() + ".icm"

    for pdbPath in pdbPaths:
        pdbName = basename(pdbPath).replace(".pdb", "")
        # This is to reflect the fact that ICM turns hyphens to underscores
        # when a file is loaded in
        pdbName = pdbName.replace("-", "_")

        # Copy the superimpose ICM script to a temp location
        shutil.copy(script, tempScript)

        # Make the sed changes
        os.system("sed -e 's|PDB_PATH_1|" + templatePath + "|g' " +
                  tempScript + " -i")
        os.system("sed -e 's|PDB_PATH_2|" + pdbPath + "|g' " +
                  tempScript + " -i")
        os.system("sed -e 's|PDB_NAME_1|" + templateName + "|g' " +
                  tempScript + " -i")
        os.system("sed -e 's|PDB_NAME_2|" + pdbName + "|g' " +
                  tempScript + " -i")

        print pdbName
        # Execute the temp script
        try:
            check_output(icm + " -s " + tempScript, stderr=STDOUT, shell=True)
        except CalledProcessError, e:
            print "\n Error during the ICM script execution"
            print e.output
            sys.exit()

        # Remove the temp script
        os.remove(tempScript)


if __name__ == "__main__":
    main()
