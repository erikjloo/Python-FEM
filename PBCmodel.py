# Import Standard Libraries
import scipy as np
import matplotlib
import matplotlib.pyplot as plt
from tkinter import Tk, filedialog

#Import Local Libraries
from solidModel import SolidModel
from shapes import Line2, Tri3

class PBCmodel(SolidModel):
            
    # Public:
    def __init__(self, path, rank=2):
        self.coords = []
        self.inod = -1
        self.nnod = 0

        self.connectivity = []
        self.iele = -1
        self.nele = 0

        self.props = []
        self.Phys = {}
        self.nPhys = 0

        self.types = []
        self.idof = 0
        self.ndof = 0

        # Read Mesh
        self.readMesh(path)
        self.rank = rank

        # Get dofs
        self.dofspace = np.empty((self.nnod, rank))
        self.dofspace[:] = np.nan
        self.addTypes(["u", "v"])
        self.shape = Tri3()

        # Add dofs
        self.addDofs(range(self.nnod), self.types)
        # Add displacement doftypes

        # Create boundary element
        self.bshape = Line2()
        self.nIP_ = self.bshape.nIP
        self.nnod_ = self.bshape.nnod
        self.localrank_ = self.bshape.ndim

        self.__boundingBox()
        self.__setTolerances()
        self.__findBndNodes()
        self.__sortBndNodes()
        self.__createTractionMesh()

    def getBndNodes(self):
        """ Output: bndNodes = [ xmin, xmax, ymin, ymax ] """
        return self.bndNodes

    def getTrNodes(self):
        """ Output: trNodes = [ xmin, ymin ] """
        return self.trNodes

    def plotBoundary(self):
        """ Plots the boundary nodes """
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111)

        for face in self.bndNodes:
            coords = self.getCoords(face)
            ax.plot(coords[:, 0], coords[:, 1],
                    marker='o', linewidth=0.3, markersize=3, color='k')
            ax.plot([coords[3, 0]], [coords[3, 1]],
                    marker='o', markersize=6, color="green")
            ax.plot([coords[7, 0]], [coords[7, 1]],
                    marker='o', markersize=6, color="red")
        plt.show()

    def __boundingBox(self):
        """ Sets box = [xmin, xmax, ymin, ymax] """
        self.box = np.ones(4)
        self.box[0:2] *= self.coords[0][0]
        self.box[2:4] *= self.coords[0][1]
        for coord in self.coords:
            if coord[0] > self.box[1]:
                self.box[1] = coord[0]
            if coord[0] < self.box[0]:
                self.box[0] = coord[0]
            if coord[1] > self.box[3]:
                self.box[3] = coord[1]
            if coord[1] < self.box[2]:
                self.box[2] = coord[1]

    def __setTolerances(self):
        """ Sets tolerances for finding boundary nodes """
        self.xtol = abs(self.box[1] - self.box[0])/1000000
        self.ytol = abs(self.box[3] - self.box[2])/1000000

    def __findBndNodes(self):
        """ Finds boundary nodes according to the set tolerances """
        self.bndNodes = [[], [], [], []]
        for inod, coord in enumerate(self.coords):
            if abs(coord[0] - self.box[0]) < self.xtol:
                self.bndNodes[0].append(inod)
            if abs(coord[0] - self.box[1]) < self.xtol:
                self.bndNodes[1].append(inod)
            if abs(coord[1] - self.box[2]) < self.ytol:
                self.bndNodes[2].append(inod)
            if abs(coord[1] - self.box[3]) < self.ytol:
                self.bndNodes[3].append(inod)

    def __sortBndNodes(self):
        """ Sorts boundary nodes in ascending x or y coordinate """
        for face, bndFace in enumerate(self.bndNodes):
            bndFace = self.__sortBndFace(bndFace)
            self.bndNodes[face] = bndFace
        
    def __sortBndFace(self, bndFace):
        """ Sorts boundary face in ascending x or y coordinate """
        for i in range(len(bndFace)):
            for j in range(len(bndFace) - 1 - i):
                c1 = self.coords[bndFace[j+1]]
                c0 = self.coords[bndFace[j]]
                if c0[1] > c1[1] or c0[0] > c1[0]:
                    bndFace[j + 1], bndFace[j] = bndFace[j], bndFace[j+1]
        return bndFace

    def __createTractionMesh(self):

        self.trNodes = [[], []]
        # loop over face pairs (ix)
        for ix in range(self.rank):
            inodes = self.bndNodes[ix*2]
            jnodes = self.bndNodes[ix*2 + 1]

            # loop over node indices in inodes
            for inod in inodes:
                coord = self.getCoords(inod)
                knod = self.addNode(coord)
                self.trNodes[ix].append(knod)
            
            # loop over node indices in jnodes
            for jnod in jnodes:
                coord = self.getCoords(jnod)
                knod = self.addNode(coord)
                self.trNodes[ix].append(knod)
            
            self.trNodes[ix] = self.__sortBndFace(self.trNodes[ix])


            

    def getCornerNodes(self):
        pass
    def augmentKsys(self):
        pass


#===========================================================================
#   Example
#===========================================================================


if __name__ == '__main__':

    mesh = PBCmodel("rve.msh",rank=2)
    bndNodes = mesh.getBndNodes()
    trNodes = mesh.getTrNodes()
    print("\n Boundary Nodes: \n", bndNodes)
    print("\n Traction Nodes: \n", trNodes)
    mesh.plotBoundary()

