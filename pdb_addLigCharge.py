#!/usr/bin/env python

# https://github.com/thomas-coudrat/toolbx_pdb
# Thomas Coudrat <thomas.coudrat@gmail.com>

import argparse
import glob
import os
import shutil
import sys
from subprocess import check_output, STDOUT, CalledProcessError
import time

def main():
    """
    Executing script
    """
    # Parsing arguemnts
    pdb_dir = parsing()

    # Generate global paths for ICM executable and ICM script
    getGlobalPaths()

    # Get the .PDB paths from the pdb_dir
    pdb_paths = sorted(glob.glob(pdb_dir + "/*.pdb"))
    print("\n###")
    print("### Adding ligand charges to PDBs in directory: " + pdb_dir)
    print("###")

    # Write a log of what was process in the corresponding directory
    log(pdb_paths, pdb_dir)

    # Run the process on PDB directory
    process_dir(pdb_paths, pdb_dir)

def process_dir(pdb_paths, pdb_dir):
    """
    Go through directory and run command on each PDB file
    """
    if len(pdb_paths) > 0:
        print("\nOpening, adding ligand charge and overwriting:")
        for pdb_path in pdb_paths:
            add_charge_and_overwrite(pdb_path)
    else:
        print("\nNo PDB files found in ", pdb_dir)

def add_charge_and_overwrite(pdb_path):
    """
    Read in .ICM script, modify for current pdb_file, execute
    """
    # Get the PDB name as it will be created by ICM
    pdb_name = os.path.basename(pdb_path).replace("-", "_").replace(".pdb", "")
    # Create a temporary script
    temp_script = os.path.basename(pdb_path).replace(".pdb", ".icm")

    # Make a temporary copy of the script to be modified
    shutil.copy(icm_template_script, temp_script)

    # Add pdb_path and pdb_name to the temporary ICM script
    path_sed = "sed -e 's|PDB_PATH|" + pdb_path + "|g' " + temp_script + " -i"
    #print(path_sed, "\n")
    os.system(path_sed)
    name_sed = "sed -e 's|PDB_NAME|" + pdb_name + "|g' " + temp_script + " -i"
    #print(name_sed, "\n")
    os.system(name_sed)

    # Print conformation name to terminal
    print("\t", pdb_path)

    # Execute the temp script
    try:
        command = icm_executable + " -s " + temp_script
        #print(command)
        check_output(command, stderr=STDOUT, shell=True)
    except CalledProcessError as e:
        print("\n Error during the ICM script execution")
        print(e.output)
        sys.exit()

    # Remove the temp script
    os.remove(temp_script)

def getGlobalPaths():
    """
    Figure the paths to the ICM executable and the ICM charge script. Generate
    global paths for them.
    """

    # Get the directory path where this Python script was executed
    pathExec = os.path.dirname(os.path.realpath(__file__))
    # The ICM script is located in that directory
    global icm_template_script
    icm_template_script = pathExec + "/ligCharge.icm"

    # Get environment variable. It returns None if not defined on the system.
    icmHome = os.environ.get('ICMHOME')

    # Return path to executable if the environment variable was found
    if icmHome == None:
        "The ICMHOME environment variable must be set for your system. Exiting."
        sys.exit()
    else:
        global icm_executable
        icm_executable = icmHome + "/icm64"

def log(pdb_paths, pdb_dir):
    """
    Write a log text file that contains the following information: time/date,
    template and list of conformations processed.
    """
    # Getting time information
    date_str = time.strftime("%H:%M:%S")
    time_str = time.strftime("%d/%m/%Y")
    # Writing log file in pdb_dir
    with open(pdb_dir + "/ligCharge.log", "w") as log_file:
        log_file.write("###\n")
        log_file.write("### Log of conformations processed for lig charges\n")
        log_file.write("###\n")
        log_file.write("Date & Time: " + date_str + " -- " + time_str + "\n\n")
        log_file.write("Processing " + str(len(pdb_paths)) + \
                       " conformations from directory: " + pdb_dir + "\n")
        for pdb in pdb_paths:
            log_file.write("\t" + os.path.basename(pdb) + "\n")

def parsing():
    """
    Parsing arguments
    """
    descr = "Add ligand charge to PDB files in directory, overwrite"
    descr_pdbDir = "Directory containing PDB files"
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument("pdb_dir", help=descr_pdbDir)

    args = parser.parse_args()

    return args.pdb_dir

if __name__ == "__main__":
    main()
