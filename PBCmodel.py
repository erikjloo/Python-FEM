# Import Standard Libraries
import scipy as np
import matplotlib.pyplot as plt
from copy import deepcopy
from tkinter import Tk, filedialog

#Import Local Libraries
from solidModel import SolidModel
from properties import Properties
from shapes import Line2, Tri3


#===========================================================================
#   PBCmodel
#===========================================================================


class PBCmodel(SolidModel):
    """ Periodic Boundary Conditions Model

    Static Members:
        __type__ = "Input is not list or array!"
        __type_int__ = "Input is not int!"
        __type_str__ = "Input is not str!"
        __renumber__ = "Erasing dofs: Dof numbers will be renumbered!"
        __rank_error__ = "Rank has to be 1, 2 or 3!"
        
    Instance Members:
        coords = list of nodal coordinates
        inod = last node index = nnod - 1
        nnod = number of nodes

        connectivity = list of element connectivities
        iele = last element index = nele - 1
        nele = number of elements

        groups = list of elements in each physical group
        groupNames = names of physical groups
        ngroups = number of physical groups

        dofspace = array of dof indices (idofs)
        types = list of dof type names
        idof = last dof index = ndof - 1
        ndof = number of dofs

        rank = number of dimensions
        U_doftypes = displacement dof types
        T_doftypes = traction dof types
        shape = element shape
        bshape = boundary element shape
        nIP_ = number of integration points of boundary element
        nnod_ = number of nodes of boundary element
        localrank_ = local rank of boundary element

        bndNodes = boundary nodes = [xmin, xmax, ymin, ymax]
        trNodes = traction nodes = [xmin, ymin]
        corner0 = corner at intersection of xmin & ymin
        cornerx = corner at intersection of xmax & ymin
        cornery = corner at intersection of xmin & ymax

    Public Methods:
        PBCmodel(path, rank=2)

    Private Methods:
        __boundingBox()
        __setTolerances()
        __findBndNodes()
        __sortBndNodes()
        bndFace = __sortBndFace(bndFace, index)
        __findCornerNodes()
        __createTractionMesh(factor)
        trFace = __coarsenMesh(trFace, index, factor)
        __getTractionMeshNodes
        __augmentMatrix
    """

    # Public:
    def __init__(self, props, name, rank=2):

        # Call the SolidModel constructor
        SolidModel.__init__(self, props, name, rank)
        self.__initialize()
        
    #-----------------------------------------------------------------------
    #   __initialize
    #-----------------------------------------------------------------------

    def __initialize(self):
        
        # Add displacement dofs
        self.U_doftypes = ["u", "v"]
        self.addTypes(self.U_doftypes)
        self.addDofs(range(self.nnod), self.U_doftypes)

        # Add traction dofs
        self.T_doftypes = ["tx","ty"]
        self.addTypes(self.T_doftypes)

        # Get specimen dimensions
        print("\n Box dimensions: \n")
        self.__boundingBox()
        self.__setTolerances()

        # Get boundary nodes
        print("\n Boundary Nodes: \n")
        self.__findBndNodes()
        self.__sortBndNodes()

        # Find corner nodes
        print("\n Corner Nodes: \n")
        self.__findCornerNodes()

        # Create traction mesh
        print("\n Traction Nodes: \n")
        self.__createTractionMesh()

        # Create boundary element
        self.shape = Tri3()
        self.bshape = Line2()
        self.nIP_ = self.bshape.nIP
        self.nnod_ = self.bshape.nnod
        self.localrank_ = self.bshape.ndim

    #-----------------------------------------------------------------------
    #   __plotBoundary
    #-----------------------------------------------------------------------

    def plotBoundary(self):
        """ Plots the boundary nodes """

        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111)

        for bndFace in self.bndNodes:
            coords = self.getCoords(bndFace)
            ax.plot(coords[:, 0], coords[:, 1],
                    marker='o', linewidth=0.3, markersize=3, color='k')
            ax.plot([coords[3, 0]], [coords[3, 1]],
                    marker='o', markersize=6, color="green")
            ax.plot([coords[7, 0]], [coords[7, 1]],
                    marker='o', markersize=6, color="red")
        
        #self.trNodes[0] = [ 2361, 2363, 2366, 2369, 2372, 2374, 2377, 2380, 2383, 2386, 2389, 2392, 2395, 2398, 2401, 2448 ]
        #self.trNodes[1] = [ 2449, 2450, 2452, 2455, 2458, 2460, 2462, 2465, 2468, 2470, 2472, 2475, 2477, 2479, 2520 ]
        
        for ix, trFace in enumerate(self.trNodes):
            coords = self.getCoords(trFace)
            coords[:,ix] -= self.dx[ix]/10
            ax.plot(coords[:,0], coords[:,1],
                    marker='o', linewidth=0.3, markersize=3, color='blue')
        
        plt.show()

    #-----------------------------------------------------------------------
    #   __boundingBox
    #-----------------------------------------------------------------------

    def __boundingBox(self):
        """ Sets box = [xmin, xmax, ymin, ymax] """
        
        # Create space for box
        self.box = np.ones(2*self.rank)
        self.box[0:2] *= self.coords[0][0]
        self.box[2:4] *= self.coords[0][1]

        # Create space for dx
        self.dx = np.zeros(self.rank)
        
        # Find specimen coordinates
        for coord in self.coords:
            for ix in range(self.rank):
                if coord[ix] < self.box[2*ix]:
                    self.box[2*ix] = coord[ix]
                if coord[ix] > self.box[2*ix+1]:
                    self.box[2*ix+1] = coord[ix]
        
        # Find specimen dimensions
        for ix in range(self.rank):
            self.dx[ix] = self.box[2*ix + 1] - self.box[2*ix]
        
        # Print to verify
        print("box = ",self.box)
        print("dx = ", self.dx)

    #-----------------------------------------------------------------------
    #   __setTolerances
    #-----------------------------------------------------------------------

    def __setTolerances(self):
        """ Sets tolerances for finding boundary nodes """
        self.xtol = []
        for ix in range(self.rank):
            self.xtol.append(abs(self.box[2*ix + 1] - self.box[2*ix])/1000000)

    #-----------------------------------------------------------------------
    #   __findBndNodes
    #-----------------------------------------------------------------------
    
    def __findBndNodes(self):
        """ Finds bndNodes according to the set tolerances """

        # Creates space for bndNodes
        self.bndNodes = [[] for face in range(2*self.rank)]

        # Loop over every node
        for inod, coord in enumerate(self.coords):
            for ix in range(self.rank):
                if abs(coord[ix] - self.box[2*ix]) < self.xtol[ix]:
                    self.bndNodes[2*ix].append(inod)
                if abs(coord[ix] - self.box[2*ix + 1]) < self.xtol[ix]:
                    self.bndNodes[2*ix + 1].append(inod)

    #-----------------------------------------------------------------------
    #   __sortBndNodes
    #-----------------------------------------------------------------------

    def __sortBndNodes(self):
        """ Sorts bndNodes in ascending x or y coordinate """

        # Loop over faces of bndNodes
        for face, bndFace in enumerate(self.bndNodes):

            # Map face onto ix
            ix = np.floor(face/2)
            # Get correct index
            index = 1 if (ix == 0) else 0
            # Perform bubblesort on bndFace
            self.bndNodes[face] = self.__sortBndFace(bndFace, index)
            # Print to verify
            print(" bndNodes[", face, "] =", self.bndNodes[face])

    #-----------------------------------------------------------------------
    #   __sortBndFace
    #-----------------------------------------------------------------------

    def __sortBndFace(self, bndFace, index):
        """ Sorts bndFace in ascending x or y coordinate """

        # Bubblesort algorithm
        for inod in range(len(bndFace)):
            for jnod in range(len(bndFace) - 1 - inod):

                # Get nodal coordinates
                c0 = self.coords[bndFace[jnod]]
                c1 = self.coords[bndFace[jnod+1]]

                # Swap indices if necessary
                if c0[index] > c1[index]:
                    bndFace[jnod + 1], bndFace[jnod] = bndFace[jnod], bndFace[jnod + 1]

        return bndFace

    #-----------------------------------------------------------------------
    #   __findCornerNodes
    #-----------------------------------------------------------------------

    def __findCornerNodes(self):
        """ Finds the intersection of each bndFace """

        if self.rank == 2:
            self.corner0 = self.bndNodes[0][0]
            self.cornerx = self.bndNodes[1][0]
            self.cornery = self.bndNodes[3][0]
            print(" corner0 = ", self.corner0)
            print(" cornerx = ", self.cornerx)
            print(" cornery = ", self.cornery)
        elif self.rank == 3:
            raise NotImplementedError(" Not yet implemented. ")

    #-----------------------------------------------------------------------
    #   __createTractionMesh
    #-----------------------------------------------------------------------

    def __createTractionMesh(self):
        """ Maps all nodes onto xmin and ymin to create traction mesh """
        
        #---------------------------------------------------------------------------
        # Part 1: Find the smallest element dimensions along the x and y coordinates
        #---------------------------------------------------------------------------
        
        # smallest element dimensions
        self.dx0 = deepcopy(self.dx) 

        # Loop over faces of bndNodes
        for face, bndFace in enumerate(self.bndNodes):
            
            # Map face onto ix
            ix = np.floor(face/2)

            # Get correct index
            index = 1 if (ix == 0) else 0
            
            # Loop over indices of bndFace
            for inod in range(len(bndFace)-1):

                # Get nodal coordinates
                c0 = self.coords[bndFace[inod]]
                c1 = self.coords[bndFace[inod+1]]

                # Calculate dx and compare to dx0
                dx = c1[index] - c0[index]
                self.dx0[index] = dx if (dx < self.dx0[index]) else self.dx0[index]

        #---------------------------------------------------------------------------
        # Part 2: Map all nodes onto xmin and ymin, reorder them and coarsen the mesh
        #---------------------------------------------------------------------------
        
        # Create space for trNodes
        self.trNodes = [[] for ix in range(self.rank)]

        # Loop over faces of trNodes
        for ix, trFace in enumerate(self.trNodes):
            inodes = self.bndNodes[2*ix]        # bndFace_min
            jnodes = self.bndNodes[2*ix + 1]    # bndFace_max

            # Loop over indices of inodes
            for inod in inodes:
                coord = self.getCoords(inod)
                trFace.append(self.addNode(coord))

            # Loop over indices of jnodes
            for inod in jnodes:
                coord = self.getCoords(inod)
                coord[ix] = self.box[2*ix]
                trFace.append(self.addNode(coord))

            # Get correct index
            index = 1 if (ix == 0) else 0

            # Sorting trFace (works only for 2D)
            trFace = self.__sortBndFace(trFace, index)

            # Coarsen trFace (works only for 2D)
            trFace = self.__coarsenMesh(trFace)

            # Assign trFace to trNodes[ix]
            self.trNodes[ix] = trFace

            # Print to verify
            print(" trNodes[", ix, "] =", self.trNodes[ix])

        # Loop over faces of trNodes
        # for trFace in self.trNodes:
            # self.addDofs(trFace,self.T_doftypes)

    #-----------------------------------------------------------------------
    #   __coarsenMesh
    #-----------------------------------------------------------------------

    def __coarsenMesh(self, trFace):
        """ Coarsens the traction mesh """
        factor = 0.3
        dx = (self.dx0[0]+self.dx0[1])/(2*factor)
        cn = self.getCoords(trFace[-1])
        
        # Loop over indices of trFace:
        for inod in range(len(trFace)):

            # Get nodal coordinates
            c0 = self.getCoords(trFace[inod])
            c1 = self.getCoords(trFace[inod+1])

            # Delete indices until c1 - c0 > dx
            while np.linalg.norm(c1 - c0) < dx:
                # Delete current index
                del trFace[inod+1]
                # Assign next node coords to c1
                c1 = self.getCoords(trFace[inod+1])

            # Check distance to last node
            if np.linalg.norm(cn - c1) < dx:
                # Delete all nodes up to but not including the last one
                del trFace[inod+1:-1]
                break

        return trFace
        
    #-----------------------------------------------------------------------
    #   __getTractionMeshNodes
    #-----------------------------------------------------------------------

    def __getTractionMeshNodes(self, x, face):
        """ Gets nodes of traction mesh element at given global x """
        
        # Map face onto ix
        ix = np.floor(face/2)

        # Implementation for two dimensions
        if self.rank == 2:
            
            # Assign trNodes[ix] onto trFace
            trFace = self.trNodes[ix]

            # Loop over indices of trFace
            for inod in range(len(trFace)-1):

                # Get coords of nodes in trFace[inod: inod+2]
                connect = trFace[inod:inod+2]
                coords = self.getCoords(connect)

                # Get correct index
                index = 1 if (ix == 0) else 0

                # Check if c0[index] < x[index] < c1[index]
                if coords[0,index] < x[index] < coords[1,index]:
                    print(coords[0,index]," < ",x[index]," < ",coords[1,index])
                    return connect

        elif self.rank == 3:
            raise NotImplementedError(" Not yet implemented. ")

    #-----------------------------------------------------------------------
    #   __augmentMatrix
    #-----------------------------------------------------------------------

    def __augmentMatrix(self):

        print("\n Augmenting Matrix: \n")
        
        # Matrix to be assembled: K[idofs, jdofs] += w[ip]*N[ip]*H[ip]
        N = np.zeros((self.rank, self.nnod_))
        H = np.zeros((self.rank, self.nnod_))
        Ke = np.zeros((self.nnod_*self.rank, self.nnod_*self.rank))

        # Loop over faces of bndNodes
        for face, bndFace in enumerate(self.bndNodes):

            # Loop over indices of bndFace
            for inod in range(len(bndFace)-1):

                # Get idofs, w and N from displacement mesh
                connect = bndFace[inod:inod+2]
                idofs = self.getDofIndices(connect, ['u', 'v'])
                coords = self.getCoords(connect)
                w = self.bshape.getIntegrationWeights(coords)
                n = self.bshape.getShapeFunctions()
                X = self.bshape.getGlobalIntegrationPoints(coords)

                # Get jdofs from traction mesh
                connect = self.__getTractionMeshNodes(X[0], face)
                jdofs = self.getDofIndices(connect, ['tx', 'ty'])
                coords = self.getCoords(connect)
                
                for ip in range(self.nIP_):
                    
                    # Assemble N matrix
                    N[0,0] = N[1,1] = n[ip,0]
                    N[0,2] = N[1,2] = n[ip,1]

                    # Assemble H matrix
                    xi = self.bshape.getLocalPoint(X[ip], coords)
                    h = self.bshape.evalShapeFunctions(xi)
                    H[0,0] = H[1,1] = h[0]
                    H[0,2] = H[1,2] = h[1]

                    Ke += w[ip] * H.transpose() @ N
                    KeT += w[ip] * N.transpose() @ H




#===========================================================================
#   Example
#===========================================================================


if __name__ == '__main__':

    file = "Examples/square.pro"
    props = Properties()
    props.parseFile(file)

    model = PBCmodel(props, "model.matrix", rank=2)
    #model.initialize()


    model.plotBoundary()
