# Import Standard Libraries
import scipy as np
import matplotlib.pyplot as plt
from copy import deepcopy

# Import Local Libraries
from models import Model
from shapes import Line2


#===========================================================================
#   PBCmodel
#===========================================================================


class PBCmodel(Model):
    """ Periodic Boundary Conditions Model
        
    Instance Members:
        rank_ = number of dimensions
        U_doftypes_ = displacement dof types
        T_doftypes_ = traction dof types
        bshape_ = boundary element shape
        nIP_ = number of integration points of boundary element
        nnod_ = number of nodes of boundary element
        localrank_ = local rank of boundary element

        bndNodes = boundary nodes = [xmin, xmax, ymin, ymax]
        trNodes = traction nodes = [xmin, ymin]
        corner0 = corner at intersection of xmin & ymin
        cornerx = corner at intersection of xmax & ymin
        cornery = corner at intersection of xmin & ymax

    Public Methods:
        PBCmodel(name, props, mesh)
        get_Matrix_0(mesh, mbuild, fint)
        get_Ext_Vector(fext)
        get_Constraints(cons, mesh)

    Private Methods:
        __boundingBox(mesh)
        __setTolerances()
        __findBndNodes(mesh)
        __sortBndNodes(mesh)
        __sortBndFace(mesh, bndFace, index)
        __findCornerNodes()
        __createTractionMesh(mesh)
        __coarsenMesh(mesh, trFace)
        __getTractionMeshNodes
    """

    # Public:

    #-----------------------------------------------------------------------
    #   initialize
    #-----------------------------------------------------------------------

    def __init__(self, name, props, mesh):

        # Model name
        self.name = name
        self.rank_ = mesh.rank
        self.factor = props.get("coarsenFactor")
        self.strain = props.get("imposedStrain")

        # Add displacement doftypes
        self.U_doftypes_ = ["u", "v"]
        mesh.addTypes(self.U_doftypes_)

        # Add traction doftypes
        self.T_doftypes_ = ["tx", "ty"]
        mesh.addTypes(self.T_doftypes_)

        # Get specimen dimensions
        print("\n Box dimensions: \n")
        self.__boundingBox(mesh)
        self.__setTolerances()

        # Get boundary nodes
        print("\n Boundary Nodes: \n")
        self.__findBndNodes(mesh)
        self.__sortBndNodes(mesh)

        # Find corner nodes
        print("\n Corner Nodes: \n")
        self.__findCornerNodes()

        # Create traction mesh
        print("\n Traction Nodes: \n")
        self.__createTractionMesh(mesh)

        # Create boundary element
        self.bshape_ = Line2()
        self.nIP_ = self.bshape_.nIP
        self.nnod_ = self.bshape_.nnod
        self.localrank_ = self.bshape_.ndim

    #-----------------------------------------------------------------------
    #   plotBoundary
    #-----------------------------------------------------------------------

    def plotBoundary(self, mesh):
        """ Plots the boundary nodes """

        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111)

        for bndFace in self.bndNodes:
            coords = mesh.getCoords(bndFace)
            ax.plot(coords[:, 0], coords[:, 1],
                    marker='o', linewidth=0.3, markersize=3, color='k')
            ax.plot([coords[3, 0]], [coords[3, 1]],
                    marker='o', markersize=6, color="green")
            ax.plot([coords[7, 0]], [coords[7, 1]],
                    marker='o', markersize=6, color="red")

        for ix, trFace in enumerate(self.trNodes):
            coords = mesh.getCoords(trFace)
            coords[:, ix] -= self.dx[ix]/10
            ax.plot(coords[:, 0], coords[:, 1],
                    marker='o', linewidth=0.3, markersize=3, color='blue')

        plt.show()

    #-----------------------------------------------------------------------
    #   __boundingBox
    #-----------------------------------------------------------------------

    def __boundingBox(self, mesh):
        """ Sets box = [xmin, xmax, ymin, ymax] """

        # Create space for box
        self.box = np.ones(2*self.rank_)
        self.box[0:2] *= mesh.coords[0][0]
        self.box[2:4] *= mesh.coords[0][1]

        # Create space for dx
        self.dx = np.zeros(self.rank_)

        # Find specimen coordinates
        for coord in mesh.coords:
            for ix in range(self.rank_):
                if coord[ix] < self.box[2*ix]:
                    self.box[2*ix] = coord[ix]
                if coord[ix] > self.box[2*ix+1]:
                    self.box[2*ix+1] = coord[ix]

        # Find specimen dimensions
        for ix in range(self.rank_):
            self.dx[ix] = self.box[2*ix + 1] - self.box[2*ix]

        # Print to verify
        print("box = ", self.box)
        print("dx = ", self.dx)

    #-----------------------------------------------------------------------
    #   __setTolerances
    #-----------------------------------------------------------------------

    def __setTolerances(self):
        """ Sets tolerances for finding boundary nodes """
        self.xtol = []
        for ix in range(self.rank_):
            self.xtol.append(abs(self.box[2*ix + 1] - self.box[2*ix])/1000000)

    #-----------------------------------------------------------------------
    #   __findBndNodes
    #-----------------------------------------------------------------------

    def __findBndNodes(self, mesh):
        """ Finds bndNodes according to the set tolerances """

        # Creates space for bndNodes
        self.bndNodes = [[] for face in range(2*self.rank_)]

        # Loop over every node
        for inod, coord in enumerate(mesh.coords):
            for ix in range(self.rank_):
                if abs(coord[ix] - self.box[2*ix]) < self.xtol[ix]:
                    self.bndNodes[2*ix].append(inod)
                if abs(coord[ix] - self.box[2*ix + 1]) < self.xtol[ix]:
                    self.bndNodes[2*ix + 1].append(inod)

    #-----------------------------------------------------------------------
    #   __sortBndNodes
    #-----------------------------------------------------------------------

    def __sortBndNodes(self, mesh):
        """ Sorts bndNodes in ascending x or y coordinate """

        # Loop over faces of bndNodes
        for face in range(len(self.bndNodes)):

            # Map face onto ix
            ix = int(np.floor(face/2))
            # Get correct index
            index = 1 if (ix == 0) else 0
            # Perform bubblesort on bndFace
            self.__sortBndFace(mesh, self.bndNodes[face], index)
            # Print to verify
            print(" bndNodes[", face, "] =", self.bndNodes[face])

    #-----------------------------------------------------------------------
    #   __sortBndFace
    #-----------------------------------------------------------------------

    def __sortBndFace(self, mesh, bndFace, index):
        """ Sorts bndFace in ascending x or y coordinate """

        # Bubblesort algorithm
        for inod in range(len(bndFace)):
            for jnod in range(len(bndFace) - 1 - inod):

                # Get nodal coordinates
                c0 = mesh.coords[bndFace[jnod]]
                c1 = mesh.coords[bndFace[jnod+1]]

                # Swap indices if necessary
                if c0[index] > c1[index]:
                    bndFace[jnod + 1], bndFace[jnod] = bndFace[jnod], bndFace[jnod + 1]

    #-----------------------------------------------------------------------
    #   __findCornerNodes
    #-----------------------------------------------------------------------

    def __findCornerNodes(self):
        """ Finds the intersection of each bndFace """

        self.corner = []
        self.corner0 = self.bndNodes[0][0]
        if self.rank_ == 2:
            self.corner.append(self.bndNodes[1][0])
            self.corner.append(self.bndNodes[3][0])
            print(" corner0 = ", self.corner0)
            print(" cornerx = ", self.corner[0])
            print(" cornery = ", self.corner[1])
        elif self.rank_ == 3:
            raise NotImplementedError(" Not yet implemented. ")

    #-----------------------------------------------------------------------
    #   __createTractionMesh
    #-----------------------------------------------------------------------

    def __createTractionMesh(self, mesh):
        """ Maps all nodes onto xmin and ymin to create traction mesh """

        #---------------------------------------------------------------------------
        # Part 1: Find the smallest element dimensions along the x and y coordinates
        #---------------------------------------------------------------------------

        # smallest element dimensions
        self.dx0 = deepcopy(self.dx)

        # Loop over faces of bndNodes
        for face, bndFace in enumerate(self.bndNodes):

            # Map face onto ix
            ix = int(np.floor(face/2))

            # Get correct index
            index = 1 if (ix == 0) else 0

            # Loop over indices of bndFace
            for inod in range(len(bndFace)-1):

                # Get nodal coordinates
                c0 = mesh.coords[bndFace[inod]]
                c1 = mesh.coords[bndFace[inod+1]]

                # Calculate dx and compare to dx0
                dx = c1[index] - c0[index]
                self.dx0[index] = dx if (
                    dx < self.dx0[index]) else self.dx0[index]

        #---------------------------------------------------------------------------
        # Part 2: Map all nodes onto xmin and ymin, reorder them and coarsen the mesh
        #---------------------------------------------------------------------------

        # Create space for trNodes
        self.trNodes = [[] for ix in range(self.rank_)]

        # Loop over faces of trNodes
        for ix, trFace in enumerate(self.trNodes):
            inodes = self.bndNodes[2*ix]        # bndFace_min
            jnodes = self.bndNodes[2*ix + 1]    # bndFace_max
            mesh.addRows(len(inodes)+len(jnodes))

            # Loop over indices of inodes
            for inod in inodes:
                coord = mesh.getCoords(inod)
                trFace.append(mesh.addNode(coord))

            # Loop over indices of jnodes
            for inod in jnodes:
                coord = mesh.getCoords(inod)
                coord[ix] = self.box[2*ix]
                trFace.append(mesh.addNode(coord))

            # Get correct index
            index = 1 if (ix == 0) else 0

            # Sorting trFace (works only for 2D)
            self.__sortBndFace(mesh, trFace, index)

            # Coarsen trFace (works only for 2D)
            self.__coarsenMesh(mesh, trFace)

            # Assign trFace to trNodes[ix]
            self.trNodes[ix] = trFace

            # Print to verify
            print(" trNodes[", ix, "] =", self.trNodes[ix])

        # Loop over faces of trNodes
        for trFace in self.trNodes:
            mesh.addDofs(trFace, self.T_doftypes_)

        # Obtain ndof --> Used in mbuild.resize(ndof)
        self.ndof = mesh.dofCount()

    #-----------------------------------------------------------------------
    #   __coarsenMesh
    #-----------------------------------------------------------------------

    def __coarsenMesh(self, mesh, trFace):
        """ Coarsens the traction mesh """

        dx = (self.dx0[0]+self.dx0[1])/(2*self.factor)
        cn = mesh.getCoords(trFace[-1])

        # Loop over indices of trFace:
        for inod in range(len(trFace)):

            # Get nodal coordinates
            c0 = mesh.getCoords(trFace[inod])
            c1 = mesh.getCoords(trFace[inod+1])

            # Delete indices until c1 - c0 > dx
            while np.linalg.norm(c1 - c0) < dx:
                # Delete current index
                del trFace[inod+1]
                # Assign next node coords to c1
                c1 = mesh.getCoords(trFace[inod+1])

            # Check distance to last node
            if np.linalg.norm(cn - c1) < dx:
                # Delete all nodes up to but not including the last one
                del trFace[inod+1:-1]
                break

    #-----------------------------------------------------------------------
    #   __getTractionMeshNodes
    #-----------------------------------------------------------------------

    def __getTractionMeshNodes(self, mesh, x, face):
        """ Gets nodes of traction mesh element at given global x """

        # Map face onto ix
        ix = int(np.floor(face/2))

        # Implementation for two dimensions
        if self.rank_ == 2:

            # Assign trNodes[ix] onto trFace
            trFace = self.trNodes[ix]

            # Loop over indices of trFace
            for inod in range(len(trFace)-1):

                # Get coords of nodes in trFace[inod: inod+2]
                connect = trFace[inod:inod+2]
                coords = mesh.getCoords(connect)

                # Get correct index
                index = 1 if (ix == 0) else 0

                # Check if c0[index] < x[index] < c1[index]
                if coords[0, index] < x[index] < coords[1, index]:
                    return connect

        elif self.rank_ == 3:
            raise NotImplementedError(" Not yet implemented. ")

    #-----------------------------------------------------------------------
    #   get_Matrix_0
    #-----------------------------------------------------------------------

    def get_Matrix_0(self, mbuild, fint, disp, mesh):

        print("\n Augmenting Matrix: \n")

        # Variables related to element on U mesh
        N = np.zeros((self.rank_, self.nnod_*self.rank_))

        # Variables related to element on T mesh
        H = np.zeros((self.rank_, self.nnod_*self.rank_))

        # Loop over faces of bndNodes
        for face, bndFace in enumerate(self.bndNodes):

            # Loop over indices of bndFace
            for inod in range(len(bndFace)-1):

                # Get idofs, w and N from displacement mesh
                connect = bndFace[inod:inod+2]
                idofs = mesh.getDofIndices(connect, self.U_doftypes_)
                coords = mesh.getCoords(connect)
                w = self.bshape_.getIntegrationWeights(coords)
                n = self.bshape_.getShapeFunctions()
                X = self.bshape_.getGlobalPoints(coords)

                # Get jdofs from traction mesh
                connect = self.__getTractionMeshNodes(mesh, X[0], face)
                jdofs = mesh.getDofIndices(connect, self.T_doftypes_)
                coords = mesh.getCoords(connect)

                # Matrix to be assembled: K[idofs, jdofs]
                Ke = np.zeros((self.nnod_*self.rank_, self.nnod_*self.rank_))

                for ip in range(self.nIP_):
          
                    # Assemble N matrix
                    N[0, 0] = N[1, 1] = n[ip][0]
                    N[0, 2] = N[1, 3] = n[ip][1]

                    # Assemble H matrix
                    xi = self.bshape_.getLocalPoint(X[ip], coords)
                    h = self.bshape_.evalShapeFunctions(xi)
                    H[0, 0] = H[1, 1] = h[0]
                    H[0, 2] = H[1, 3] = h[1]

                    # Assemble Ke
                    if face == 0 or face == 2:
                        Ke -= w[ip] * N.transpose() @ H
                    else:
                        Ke += w[ip] * N.transpose() @ H

                    # Assemble fe
                    fe = Ke.dot(disp[jdofs])

                # Add fe to fint and Ke to mbuild
                fint[idofs] += fe
                mbuild.addBlock(idofs, jdofs, Ke)
                mbuild.addBlock(jdofs, idofs, Ke.transpose())

        for ix, trFace in enumerate(self.trNodes):
    
            jdofs = mesh.getDofIndices(self.corner[ix],self.T_doftypes_)
            u_corner = disp[jdofs]
            
            for inod in range(len(trFace)-1):

                # Get jdofs, w and H from traction mesh
                connect = trFace[inod:inod+2]
                jdofs = mesh.getDofIndices(connect)
                coords = mesh.getCoords(connect)
                w = self.bshape_.getIntegrationWeights(coords)
                h = self.bshape_.getShapeFunctions()

                # Vector to be assembled: fext[jdofs]
                fe = np.zeros(self.nnod_*self.rank_)

                for ip in range(self.nIP_):

                    # Assemble H matrix
                    H[0, 0] = H[1, 1] = h[ip][0]
                    H[0, 2] = H[1, 3] = h[ip][1]

                    # Assemble fe
                    fe += w[ip] * (H.transpose() @ u_corner)

                fint[jdofs] += fe

    #-----------------------------------------------------------------------
    #   get_Ext_Vector
    #-----------------------------------------------------------------------

    def get_Ext_Vector(self, fext, mesh):

        print("\n Augmenting fext: \n")

        # Variables related to element on T mesh
        H = np.zeros((self.rank_, self.nnod_*self.rank_))
        eps = np.zeros((self.rank_, self.rank_))
        eps[0, 0] = self.strain[0]
        eps[1, 1] = self.strain[1]
        eps[0, 1] = eps[1, 0] = self.strain[2]/2

        # Loop over faces of trNodes
        for ix, trFace in enumerate(self.trNodes):
            
            u_corner = eps[ix, :]*self.dx[ix]

            # Loop over indices of trFace
            for inod in range(len(trFace)-1):

                # Assemble H matrix
                connect = trFace[inod:inod+2]
                jdofs = mesh.getDofIndices(connect, self.T_doftypes_)
                coords = mesh.getCoords(connect)
                w = self.bshape_.getIntegrationWeights(coords)
                h = self.bshape_.getShapeFunctions()

                # Vector to be assembled: fext[jdofs]
                fe = np.zeros(self.nnod_*self.rank_)

                for ip in range(self.nIP_):

                    # Assemble H matrix
                    H[0, 0] = H[1, 1] = h[ip][0]
                    H[0, 2] = H[1, 3] = h[ip][1]

                    # Assemble fe
                    fe += w[ip] * (H.transpose() @ u_corner)

                fext[jdofs] += fe

    #-----------------------------------------------------------------------
    #   get_Constraints
    #-----------------------------------------------------------------------

    def get_Constraints(self, cons, mesh):
        
        # Fix corner
        idofs = mesh.getDofIndices(self.corner,self.U_doftypes_)
        cons.addConstraints(idofs)

    #-----------------------------------------------------------------------
    #   takeAction
    #-----------------------------------------------------------------------

    def takeAction(self, action, mesh):
        if action == "plot_boundary":
            self.plotBoundary(mesh)


#===========================================================================
#   Example
#===========================================================================


if __name__ == '__main__':

    from properties import Properties
    from constraints import Constraints
    from algebra import MatrixBuilder
    from modules import InputModule
    from models import ModelFactory
    from nonlin import multistep
    from mesh import Mesh

    # Initialization
    np.set_printoptions(precision=4)
    file = "Examples/rve.pro"
    props = Properties()
    props.parseFile(file)
    mesh = Mesh()

    module = InputModule("input")
    module.init(props, mesh)

    model = ModelFactory("model", props, mesh)
    # model.takeAction("plot_boundary", mesh)

    ndof = mesh.dofCount()
    cons = Constraints(ndof)
    fext = np.zeros(ndof)
    fint = np.zeros(ndof)
    disp = np.zeros(ndof)
    mbuild = MatrixBuilder(ndof)

    model.get_Ext_Vector(fext, mesh)
    model.get_Constraints(cons, mesh)
    model.get_Matrix_0(mbuild, fint, disp, mesh)
