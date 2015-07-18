#!/usr/bin/env python

# -----------------------------------------------------------------------------
# This script checks all .pdb files in a directory, and writes to a text file
# the coordinates of the C-alphas of all those .pdb files
# It also extracts information stored in the REMARK statement at the top of
# each .pdb.
# A sub-list of residues can be submitted (for example binding pocket residues)
# so that the coordinates of only those residues are extracted and saved.
#
# Thomas Coudrat, September 2013
# -----------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d


class Principal_component_analysis:

    def __init__(self, conformations, resList=None):
        """
        Create the PCA instance
        """
        self.conformations = conformations
        self.resList = resList
        self.pcaCoordsArray = None
        self.vars_data = None

    def makePCAcoords(self):
        """
        Create the PCA data on self.conformations
        """
        # This will store the coordinates of all conformations in a list
        allConfCoords = []

        # Go through the conformations in alphanumeric order
        sortedConfNames = self.conformations.keys()
        sortedConfNames.sort()

        for confName in sortedConfNames:
            # Get the dictionary, and the file path
            conformationDict = self.conformations[confName]
            confPath = conformationDict['path']

            # Get the coords from that .pdb file
            confCoords = self.getPDBcoord(confPath)
            # print confCoords
            # Add this conformation's coordinates to the master array
            allConfCoords.append(confCoords)

        # Create an array with this data, and store it in the variable
        # self.pcaCoordsArray
        self.pcaCoordsArray = np.array(allConfCoords)

    def getPDBcoord(self, pdbPath):
        """
        Get a pdbPath, read that .pdb, and store the C-alpha data for that
        conformation
        Return an array for this conformation of the C-alpha x y z coordinates
        """
        pdbFile = open(pdbPath, "r")
        pdbLines = pdbFile.readlines()
        pdbFile.close()

        # This stores the conformations for .pdb file
        confCoords = []
        for line in pdbLines:
            ll = line.split()
            if len(ll) > 1:
                # In the ATOM description line
                if ll[0] == 'ATOM':
                    # Get all the C-alphas, and make sure of the alignment
                    # (in case there is no chain name)
                    if ll[2] == 'CA':
                        # If a residue list was provided, upon creation of the
                        # PCA instance, use it to select which residue
                        # coordinates to be saved
                        if self.resList:
                            if int(line[22:26]) in self.resList:
                                xCoord = float(line[31:39].strip())
                                yCoord = float(line[39:47].strip())
                                zCoord = float(line[47:55].strip())
                                confCoords = np.concatenate(
                                    [confCoords, [xCoord]])
                                confCoords = np.concatenate(
                                    [confCoords, [yCoord]])
                                confCoords = np.concatenate(
                                    [confCoords, [zCoord]])
                        # Otherwise save all the C-alphas found in the .pdb
                        # file
                        else:
                            xCoord = float(line[31:39].strip())
                            yCoord = float(line[39:47].strip())
                            zCoord = float(line[47:55].strip())
                            confCoords = np.concatenate(
                                [confCoords, [xCoord]])
                            confCoords = np.concatenate(
                                [confCoords, [yCoord]])
                            confCoords = np.concatenate(
                                [confCoords, [zCoord]])

        return confCoords

    def makePCAvars(self, var_to_plot):
        """
        Create a self.vars_data list of variables data, to be used for the
        plotting of PCA graphs
        """
        # Initialize the self.allVariables dict, based on the variable string
        # list passed as argument
        self.vars_data = {}
        self.vars_data[var_to_plot] = []

        # Go through the conformations in alphanumeric order
        sortedConfNames = self.conformations.keys()
        sortedConfNames.sort()
        for confName in sortedConfNames:
            # Get the dictionary, and the file path
            conformationDict = self.conformations[confName]

            # Go over the sorted variable data, and for each variable dataset,
            # add the value of the current conformation
            sorted_var_names = sorted(self.vars_data.keys())
            for varName in sorted_var_names:
                confVar = conformationDict[varName]
                self.vars_data[varName].append(confVar)

    def plotPCAfig(self, var_to_plot, labels, dim):
        """
        After the coords onto which apply PCA have beend extracted and stored
        in self.pcaCoordsArray, and..
        After the variables data have been extracted and stored in
        self.vars_all_data
        This function can be called to calculate the PCA of dimension 2 or 3,
        and call the plotting function to add subplots for each variable to
        be displayed
        """
        if self.pcaCoordsArray is not None:
            print(self.pcaCoordsArray)
            pca = PCA(n_components=dim)
            X_r = pca.fit(self.pcaCoordsArray).transform(self.pcaCoordsArray)
            # print X_r

            # Percentage of variance explained for each components
            print('explained variance ratio (first ' + str(dim) +
                  ' components): %s' % str(pca.explained_variance_ratio_))
            pc_labels = pca.explained_variance_ratio_

            # print X_r[:, 0]
            # print X_r[:, 1]
            # print X_r[:, 2]

            if self.vars_data is not None:
                if var_to_plot in self.vars_data.keys():
                    # Get the variable data
                    varData = self.vars_data[var_to_plot]
                    if dim == 2:
                        fig2D = plt.figure()
                        fig2D.set_facecolor('white')
                        fig2D.canvas.set_window_title("PCA 2D")
                        self.pcaSubplot(X_r, dim, varData, var_to_plot,
                                        pc_labels, labels, fig2D, 111)
                    if dim == 3:
                        fig3D = plt.figure()
                        fig3D.set_facecolor('white')
                        fig3D.canvas.set_window_title("PCA 3D")
                        self.pcaSubplot(X_r, dim, varData, var_to_plot,
                                        pc_labels, labels, fig3D, 111)

                    plt.show()
                else:
                    print("The variable required is not in the loaded set")
            else:
                print("There is no variable data loaded")
                print("Firt run Principal_component_analysis.makePCAvars(...)")
        else:
            print("The coords onto which apply the PCA have to be extracted")
            print("First run Principal_component_analysis.makePCAvars()")

    def pcaSubplot(self, X_r, dim, varData, varName, pc_labels, labels, fig,
                   position):
        """
        Get the Principal Component Analysis data for this set of coordinates
        The value of 'dim' specifies the number of dimensions to diplay
        Then plot the PCA data
        """
        # Choose dimention
        if dim == 2:
            ax = fig.add_subplot(position)
            scat = ax.scatter(X_r[:, 0], X_r[:, 1],
                              c=varData, marker="o", lw=0)  # cmap=plt.cm.hot
            for label, x, y in zip(labels, X_r[:, 0], X_r[:, 1]):
                ax.annotate(label, xy=(x, y), fontsize=10,
                            ha='left', va='bottom')
            ax.set_xlabel("PC1 (%.3f)" % pc_labels[0])
            ax.set_ylabel("PC2 (%.3f)" % pc_labels[1])
            # plt.title('PCA-2D ' + varName)
        if dim == 3:
            ax = fig.add_subplot(position, projection='3d')
            scat = Axes3D.scatter(ax, X_r[:, 0], X_r[:, 1], X_r[:, 2],
                                  c=varData, marker="o", lw=0)
            for label, x, y, z in zip(labels, X_r[:, 0], X_r[:, 1], X_r[:, 2]):
                if label != "":
                    x2D, y2D, _ = proj3d.proj_transform(x, y, z, ax.get_proj())
                    ax.annotate(label, xy=(x2D, y2D), fontsize=10,
                                ha='left', va='bottom')
            ax.set_xlabel("PC1 (%.3f)" % pc_labels[0])
            ax.set_ylabel("PC2 (%.3f)" % pc_labels[1])
            ax.set_zlabel("PC3 (%.3f)" % pc_labels[2])
            # plt.title('PCA-3D ' + varName)
        # Plot the colorbar refering to the variable coloring the conformation
        # dots
        cb = plt.colorbar(scat)
        cb.set_label(varName)
