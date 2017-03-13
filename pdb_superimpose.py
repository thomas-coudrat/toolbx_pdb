#!/usr/bin/env python

# Work on a conformation ensemble in .pdb format and superimpose them all onto
# a user chosen .pdb file
#
# https://github.com/thomas-coudrat/toolbx_pdb
# Thomas Coudrat <thomas.coudrat@gmail.com>

import os
import sys
import glob
import shutil
import argparse
from os.path import basename
from subprocess import check_output, STDOUT, CalledProcessError
import socket
import json
import time


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

    # Print info and write log file to the directory of conformations that were
    # superimposed
    print_and_log(templatePath, pdbDir, pdbPaths)

    # Run the superimposition, using an ICM script
    superimpose(templatePath, pdbPaths, pdbDir, icm, script)

def print_and_log(templatePath, pdbDir, pdbPaths):
    """
    Write a log text file that contains the following information: time/date,
    template and list of conformations superimposed. Also print information to
    terminal.
    """
    # Print
    print("-------------------")
    print("## Superimposing ##")
    print("-------------------")
    print("\nTemplate: ")
    print("\t" + templatePath)
    print("\nSuperimposing " + str(len(pdbPaths)) + \
          " conformations from directory: " + pdbDir)

    # Log
    date_str = time.strftime("%H:%M:%S")
    time_str = time.strftime("%d/%m/%Y")
    with open(pdbDir + "/superimposed.log", "w") as log_file:
        log_file.write("#")
        log_file.write("### Log of superimposed conformations\n")
        log_file.write("#")
        log_file.write("Date & Time: " + date_str + " -- " + time_str + "\n\n")
        log_file.write("Template used:\n")
        log_file.write("\t" + templatePath + "\n\n")
        log_file.write("Conformations superimposed:\n")
        for pdb in pdbPaths:
            log_file.write("\t" + basename(pdb) + "\n")

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

    # Get the directory path where this Python script was executed
    pathExec = os.path.dirname(os.path.realpath(__file__))
    # The ICM script is located in that directory
    script = pathExec + "/super.icm"

    # Get environment variable. It returns None if not defined on the system.
    icmHome = os.environ.get('ICMHOME')

    # Return path to executable if the environment variable was found
    if icmHome == None:
        "The ICMHOME environment variable must be set for your system. Exiting."
        sys.exit()
    else:
        icm = icmHome + "/icm64"

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

        # Print conformation name to terminal
        print("\t" + pdbName)

        # Execute the temp script
        try:
            check_output(icm + " -s " + tempScript, stderr=STDOUT, shell=True)
        except CalledProcessError as e:
            print("\n Error during the ICM script execution")
            print(e.output)
            sys.exit()

        # Remove the temp script
        os.remove(tempScript)
    print("")


if __name__ == "__main__":
    main()
