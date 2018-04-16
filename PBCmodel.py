# Import Standard Libraries
import scipy as np
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

        # Create interior element
        if rank == 2:
            self.shape = Tri3()

        # Create boundary element
        if rank == 2:
            self.bshape = Line2()
        elif rank == 3:
            self.bshape = Tri3()

        self.nIP_ = self.bshape.nIP
        self.nnod_ = self.bshape.nnod
        self.localrank_ = self.bshape.ndim

        # Create dofspace
        self.dofspace = np.empty((self.nnod, rank))
        self.dofspace[:] = np.nan

        # Add displacement dofs
        self.U_doftypes = ["u", "v"]
        if rank == 3:
            self.U_doftypes.append(["w"])
        self.addTypes(self.U_doftypes)
        self.addDofs(range(self.nnod), self.U_doftypes)

        # Add traction dofs
        self.T_doftypes = ["tx","ty"]
        if rank == 3:
            self.T_doftypes.append(["tz"])
        self.addTypes(self.T_doftypes)

        # Get boundary nodes
        self.__boundingBox()
        self.__setTolerances()
        self.__findBndNodes()
        print("\n Boundary Nodes: \n")
        self.__sortBndNodes()

        # Create traction mesh
        print("\n Traction Nodes: \n")
        self.__createTractionMesh()

        # Find corner nodes
        print("\n Corner Nodes: \n")
        self.__findCornerNodes()

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
        """ Finds bndNodes according to the set tolerances """
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
        """ Sorts bndNodes in ascending x or y coordinate """

        # Loop over faces of bndNodes
        for face, bndFace in enumerate(self.bndNodes):

            # Map face onto ix
            ix = np.floor(face/2)
            # Get correct index for sortBndFace
            index = [1 if (ix == 0) else 0][0]
            # Sorting bndFace
            bndFace = self.__sortBndFace(bndFace, index)
            self.bndNodes[face] = bndFace
            # Print to verify
            print(" bndNodes[", face, "] =", self.bndNodes[face])

    def __sortBndFace(self, bndFace, index):
        """ Sorts bndFace in ascending x or y coordinate """

        # Bubblesort algorithm
        for i in range(len(bndFace)):
            for j in range(len(bndFace) - 1 - i):

                # Get nodal coordinates
                c1 = self.coords[bndFace[j+1]]
                c0 = self.coords[bndFace[j]]

                # Swap indices if necessary
                if c0[index] > c1[index]:
                    bndFace[j + 1], bndFace[j] = bndFace[j], bndFace[j+1]

        return bndFace

    def __findCornerNodes(self):
        """ Finds the intersection of each bndFace """

        if self.rank == 2:
            self.corner0 = list(
                set(self.bndNodes[0]).intersection(self.bndNodes[2]))
            self.cornerx = list(
                set(self.bndNodes[1]).intersection(self.bndNodes[2]))
            self.cornery = list(
                set(self.bndNodes[0]).intersection(self.bndNodes[3]))
            print(" corner0 = ", self.corner0)
            print(" cornerx = ", self.cornerx)
            print(" cornery = ", self.cornery)
        elif self.rank == 3:
            raise NotImplementedError(" Not yet implemented ")

    def __createTractionMesh(self):
        """ Maps all nodes onto xmin and ymin to create traction mesh """

        self.trNodes = [[], []]
        # Loop over face pairs (ix)
        for ix in range(self.rank):
            inodes = self.bndNodes[ix*2]
            jnodes = self.bndNodes[ix*2 + 1]

            # Loop over node indices in inodes
            for inod in inodes:
                coord = self.getCoords(inod)
                knod = self.addNode(coord)
                self.trNodes[ix].append(knod)

            # Loop over node indices in jnodes
            for jnod in jnodes:
                coord = self.getCoords(jnod)
                knod = self.addNode(coord)
                self.trNodes[ix].append(knod)

            # Sorting trFace (workds only for 2D)
            index = [1 if (ix == 0) else 0][0]
            self.trNodes[ix] = self.__sortBndFace(self.trNodes[ix], index)

            # Print to verify
            print(" trNodes[", ix, "] =", self.trNodes[ix])

    def __getTractionMeshNodes(self, x, face):
        """ Gets nodes of traction mesh element at given global x """
        
        # Map face index onto ix index
        ix = np.floor(face/2)

        # Implementation for two dimensions
        if self.rank == 2:
            
            # Assign trNodes[ix] onto trFace
            trFace = trNodes[ix]

            # Loop over node indices in trFace
            for inod in range(len(trFace)-1):

                # Get coords of nodes in trFace[inod: inod+2]
                connect = trFace[inod:inod+2]
                coords = self.getCoords(connect)

                # Index = 1 if ix == 0, else index = 0
                index = [1 if (ix == 0) else 0][0]
                


        pass

    def __augmentMatrix(self):
        print("\n Augmenting Matrix: \n")

        # Loop over faces of bndNodes
        for face in range(2*self.rank):

            # Assign bndNodes[face] to bndFace
            bndFace = self.bndNodes[face]

            # Loop over indices in bndFace
            for inod in range(len(bndFace)-1):

                # Get idofs, w and N from displacement mesh
                connect = bndFace[inod:inod+2]
                idofs = self.getDofIndices(connect, ['u', 'v'])
                coords = self.getCoords(connect)
                w = self.bshape.getIntegrationWeights(coords)
                N = self.bshape.getShapeFunctions()
                X = self.bshape.getGlobalIntegrationPoints(coords)

                # Get jdofs and H from traction mesh
                jnodes = __getTractionMeshNodes(X[0], face)
                jdofs = self.getDofIndices(jnodes, ['tx', 'ty'])


#===========================================================================
#   Example
#===========================================================================


if __name__ == '__main__':

    mesh = PBCmodel("Examples/rve.msh", rank=2)

    mesh.plotBoundary()
