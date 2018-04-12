# Import Standard Libraries
import scipy as np
import matplotlib
import matplotlib.pyplot as plt
from tkinter import Tk, filedialog

#Import Local Libraries
from mesh import Mesh
#sfrom itemset import NodeSet, ElementSet

class PBCmodel(Mesh):
            

    def getBoundaryNodes(self):
        """ Output: boundaryNodes = [ xmin, xmax, ymin, ymax ] """
        self.__boundingBox()
        self.__setTolerances()
        self.__findBoundaryNodes()
        self.__sortBoundaryNodes()
        return self.boundaryNodes

    def plotBoundary(self):
        """ Plots the boundary nodes """
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111)

        for row in self.boundaryNodes:
            coords = self.getCoords(row)
            ax.plot(coords[:, 0], coords[:, 1],
                    marker='o', linewidth=0.3, markersize=3, color='k')
            ax.plot([coords[3, 0]], [coords[3, 1]],
                    marker='o', markersize=6, color="green")
            ax.plot([coords[7, 0]], [coords[7, 1]],
                    marker='o', markersize=6, color="red")
        plt.show()

    def __boundingBox(self):
        """ Sets xmin, xmax, ymin and ymax coordinates """
        coord = self.coords[0]
        self.xmin = self.xmax = coord[0]
        self.ymin = self.ymax = coord[1]
        for coord in self.coords:
            if coord[0] > self.xmax:
                self.xmax = coord[0]
            if coord[0] < self.xmin:
                self.xmin = coord[0]
            if coord[1] > self.ymax:
                self.ymax = coord[1]
            if coord[1] < self.ymin:
                self.ymin = coord[1]

    def __setTolerances(self):
        """ Sets tolerances for finding boundary nodes """
        self.xtol = abs(self.xmax - self.xmin)/1000000
        self.ytol = abs(self.ymax - self.ymin)/1000000

    def __findBoundaryNodes(self):
        """ Finds boundary nodes according to the set tolerances """
        self.boundaryNodes = [[], [], [], []]
        for inod, coord in enumerate(self.coords):
            if abs(coord[0] - self.xmin) < self.xtol:
                self.boundaryNodes[0].append(inod)
            if abs(coord[0] - self.xmax) < self.xtol:
                self.boundaryNodes[1].append(inod)
            if abs(coord[1] - self.ymin) < self.ytol:
                self.boundaryNodes[2].append(inod)
            if abs(coord[1] - self.ymax) < self.ytol:
                self.boundaryNodes[3].append(inod)

    def __sortBoundaryNodes(self):
        """ Sorts boundary nodes in ascending x or y coordinate """
        for irow, row in enumerate(self.boundaryNodes):
            row = self.__sortBoundaryRow(row)
            self.boundaryNodes[irow] = row
        
    def __sortBoundaryRow(self, row):
        """ Sorts boundary row in ascending x or y coordinate """
        for i in range(len(row)):
            for j in range(len(row) - 1 - i):
                c1 = self.coords[row[j+1]]
                c0 = self.coords[row[j]]
                if c0[1] > c1[1] or c0[0] > c1[0]:
                    row[j + 1], row[j] = row[j], row[j+1]
        return row

    def getCornerNodes(self):
        pass
    def augmentKsys(self):
        pass


#===========================================================================
#   Example
#===========================================================================


if __name__ == '__main__':

    mesh = PBCmodel()
    mesh.readMesh("rve.msh")
    boundaryNodes = mesh.getBoundaryNodes()
    print(boundaryNodes)
    mesh.plotBoundary()

