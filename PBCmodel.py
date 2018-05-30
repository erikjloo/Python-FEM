# Import Standard Libraries
import scipy as np
import matplotlib.pyplot as plt
from copy import deepcopy

# Import Local Libraries
from models import Model
from shapes import ShapeFactory
from algebra import norm


#===========================================================================
#   PBCmodel
#===========================================================================


class PBCmodel(Model):
    """ Periodic Boundary Conditions Model
        
    Instance Members:
        rank = number of dimensions
        U_doftypes = displacement dof types
        T_doftypes = traction dof types
        bshape = boundary element shape
        nIP = number of integration points of boundary element
        nnod = number of nodes of boundary element
        localrank = local rank of boundary element

        bndNodes = boundary nodes = [xmin, xmax, ymin, ymax]
        trNodes = traction nodes = [xmin, ymin]
        corner0 = corner at xmin & ymin
        corner = [cornerx, cornery]

    Public Methods:
        PBCModel(name, props, mesh)
        get_Matrix_0(mbuild, fint, disp, mesh)
        get_Ext_Vector(fext, mesh)
        get_Int_Vector(fint, disp, mesh)
        get_Constraints(cons, mesh)
        takeAction(action, mesh)

    Private Methods:
        __boundingBox(mesh)
        __setTolerances()
        __findBndNodes(mesh)
        __sortBndNodes(mesh)
        __sortBndFace(mesh, bndFace, index)
        __findCornerNodes()
        __findSmallestElement(mesh)
        __createTractionMesh(mesh)
        __coarsenMesh(mesh, trFace)
        __getTractionMeshNodes(mesh, x, face)
    """

    # Public:

    #-----------------------------------------------------------------------
    #   initialize
    #-----------------------------------------------------------------------

    def __init__(self, name, props, mesh):

        # Model name
        self.name = name
        self.rank = mesh.rank
        self.factor = props.get("coarsenFactor")
        self.strain = props.get("strainRate")
        if self.factor > 1.0:
            msg = "Factor = {}. Should be less than 1".format(self.factor)
            raise ValueError(msg)

        # Add displacement doftypes
        types = ['u', 'v', 'w']
        self.U_doftypes = [types[x] for x in range(self.rank)]
        mesh.addTypes(self.U_doftypes)

        # Add traction doftypes
        types = ['tx', 'ty', 'tz']
        self.T_doftypes = [types[x] for x in range(self.rank)]
        mesh.addTypes(self.T_doftypes)

        # Get specimen dimensions
        print("    Box dimensions:")
        self.__boundingBox(mesh)
        self.__setTolerances()

        # Get boundary nodes
        print("    Boundary Nodes:")
        self.__findBndNodes(mesh)
        self.__sortBndNodes(mesh)

        # Find corner nodes
        print("    Corner Nodes:")
        self.__findCornerNodes()

        # Create traction mesh
        print("    Traction Nodes:")
        self.__findSmallestElement(mesh)
        self.__createTractionMesh(mesh)

        # Create boundary element
        self.bshape = ShapeFactory(props)
        self.localrank = self.bshape.ndim
        self.nnod = self.bshape.nnod
        self.nIP = self.bshape.nIP
        self.ndof = self.nnod*self.rank
        if self.localrank != self.rank-1:
            msg = "Shape ndim = {}. Should be {}".format(
                self.localrank, self.rank-1)
            raise ValueError(msg)

    #-----------------------------------------------------------------------
    #   __boundingBox
    #-----------------------------------------------------------------------

    def __boundingBox(self, mesh):
        """ Sets box = [xmin, xmax, ymin, ymax] """

        # Create space for box
        self.box_ = np.ones(2*self.rank)
        self.box_[0:2] *= mesh.coords[0][0]
        self.box_[2:4] *= mesh.coords[0][1]

        # Create space for dx
        self.dx_ = np.zeros(self.rank)

        # Find specimen coordinates
        for coord in mesh.coords:
            for ix in range(self.rank):
                if coord[ix] < self.box_[2*ix]:
                    self.box_[2*ix] = coord[ix]
                if coord[ix] > self.box_[2*ix+1]:
                    self.box_[2*ix+1] = coord[ix]

        # Find specimen dimensions
        for ix in range(self.rank):
            self.dx_[ix] = self.box_[2*ix + 1] - self.box_[2*ix]

        # Print to verify
        print("        box = ", self.box_)
        print("        dx = ", self.dx_)

    #-----------------------------------------------------------------------
    #   __setTolerances
    #-----------------------------------------------------------------------

    def __setTolerances(self):
        """ Sets tolerances for finding boundary nodes """
        self.xtol = []
        for ix in range(self.rank):
            self.xtol.append(
                abs(self.box_[2*ix + 1] - self.box_[2*ix])/1000000)

    #-----------------------------------------------------------------------
    #   __findBndNodes
    #-----------------------------------------------------------------------

    def __findBndNodes(self, mesh):
        """ Finds bndNodes according to the set tolerances """

        # Creates space for bndNodes
        self.bndNodes = [[] for face in range(2*self.rank)]

        # Loop over every node
        for inod, coord in enumerate(mesh.coords):
            for ix in range(self.rank):
                if abs(coord[ix] - self.box_[2*ix]) < self.xtol[ix]:
                    self.bndNodes[2*ix].append(inod)
                if abs(coord[ix] - self.box_[2*ix + 1]) < self.xtol[ix]:
                    self.bndNodes[2*ix + 1].append(inod)

    #-----------------------------------------------------------------------
    #   __sortBndNodes
    #-----------------------------------------------------------------------

    def __sortBndNodes(self, mesh):
        """ Sorts bndNodes in ascending x or y coordinate """

        # Loop over faces of bndNodes
        for face, bndFace in enumerate(self.bndNodes):

            # Map face onto ix
            ix = int(np.floor(face/2))
            # Get correct index
            index = 1 if (ix == 0) else 0
            # Perform bubblesort on bndFace
            self.__sortBndFace(mesh, bndFace, index)
            # Print to verify
            print("         bndNodes[{}] = {}".format(face, self.bndNodes[face]))

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
        if self.rank == 2:
            self.corner.append(self.bndNodes[1][0])
            self.corner.append(self.bndNodes[3][0])
            print("        corner0 = ", self.corner0)
            print("        cornerx = ", self.corner[0])
            print("        cornery = ", self.corner[1])
        elif self.rank == 3:
            raise NotImplementedError(" Not yet implemented. ")

    #-----------------------------------------------------------------------
    #   __findSmallestElement
    #-----------------------------------------------------------------------

    def __findSmallestElement(self, mesh):
        """ Finds the smallest element dimension on each bndFace """

        # smallest element dimensions
        self.dx_0 = deepcopy(self.dx_)

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
                self.dx_0[index] = dx if dx < self.dx_0[index] else self.dx_0[index]

    #-----------------------------------------------------------------------
    #   __createTractionMesh
    #-----------------------------------------------------------------------

    def __createTractionMesh(self, mesh):
        """ Maps all nodes onto xmin and ymin to create traction mesh """

        # Create space for trNodes
        self.trNodes = [[] for ix in range(self.rank)]

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
            for jnod in jnodes:
                coord = mesh.getCoords(jnod)
                coord[ix] = self.box_[2*ix]
                trFace.append(mesh.addNode(coord))

            # Get correct index
            index = 1 if (ix == 0) else 0

            # Sorting trFace (works only for 2D)
            self.__sortBndFace(mesh, trFace, index)

            # Coarsen trFace (works only for 2D)
            self.__coarsenMesh(mesh, trFace)

            # Print to verify
            print("        trNodes[{}] = {}".format(ix, self.trNodes[ix]))

        # Add dofs to traction mesh
        for trFace in self.trNodes:
            mesh.addDofs(trFace, self.T_doftypes)

    #-----------------------------------------------------------------------
    #   __coarsenMesh
    #-----------------------------------------------------------------------

    def __coarsenMesh(self, mesh, trFace):
        """ Coarsens the traction mesh """

        cn = mesh.getCoords(trFace[-1])
        dx = (self.dx_0[0]+self.dx_0[1])/(2*self.factor)

        # Loop over indices of trFace:
        for inod in range(len(trFace)):

            # Get nodal coordinates
            c0 = mesh.getCoords(trFace[inod])
            c1 = mesh.getCoords(trFace[inod+1])

            # Delete indices until c1 - c0 > dx
            while norm(c1 - c0) < dx:
                # Delete current index
                del trFace[inod+1]
                # Assign next node coords to c1
                c1 = mesh.getCoords(trFace[inod+1])

            # Check distance to last node
            if norm(cn - c1) < dx:
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
        if self.rank == 2:

            # Loop over indices of trFace
            for inod in range(len(self.trNodes[ix])-1):

                # Get coords of nodes in trFace[inod: inod+2]
                connect = self.trNodes[ix][inod:inod+2]
                coords = mesh.getCoords(connect)

                # Get correct index
                index = 1 if (ix == 0) else 0

                # Check if c0[index] < x[index] < c1[index]
                if coords[0, index] < x[index] < coords[1, index]:
                    return connect

        elif self.rank == 3:
            raise NotImplementedError(" Not yet implemented. ")

    #-----------------------------------------------------------------------
    #   get_Matrix_0
    #-----------------------------------------------------------------------

    def get_Matrix_0(self, mbuild, fint, disp, mesh):

        max_hbw = 0

        # Loop over faces of bndNodes
        for face, bndFace in enumerate(self.bndNodes):

            # Loop over indices of bndFace
            for inod in range(len(bndFace)-1):

                # Get idofs and w from U-mesh
                connect = bndFace[inod:inod+2]
                coords = mesh.getCoords(connect)
                idofs = mesh.getDofIndices(connect, self.U_doftypes)
                w = self.bshape.getIntegrationWeights(coords)
                X = self.bshape.getGlobalPoints(coords)

                # Get jdofs from T-mesh
                connect = self.__getTractionMeshNodes(mesh, X[0], face)
                jdofs = mesh.getDofIndices(connect, self.T_doftypes)
                hbw = max(jdofs) - min(idofs)
                max_hbw = hbw if hbw > max_hbw else max_hbw
                coords = mesh.getCoords(connect)

                # Matrix to be assembled: K[idofs, jdofs]
                Ke = np.zeros((self.ndof, self.ndof))

                for ip in range(self.nIP):

                    # Get N-matrix from U-mesh
                    N = self.bshape.getNmatrix(ip, ndim=self.rank)

                    # Get H-matrix from T-mesh
                    xi = self.bshape.getLocalPoint(X[ip], coords)
                    H = self.bshape.evalNmatrix(xi, ndim=self.rank)

                    # Assemble Ke
                    if face == 0 or face == 2 or face == 4:
                        Ke -= w[ip] * N.transpose() @ H
                    else:
                        Ke += w[ip] * N.transpose() @ H

                # Add Ke and KeT to mbuild
                KeT = Ke.transpose()
                if mbuild is not None:
                    mbuild.addBlock(idofs, jdofs, Ke)
                    mbuild.addBlock(jdofs, idofs, KeT)

                # Assemble U-mesh fe
                fe = Ke.dot(disp[jdofs])
                fint[idofs] += fe

                # Assemble T-mesh fe
                fe = KeT.dot(disp[idofs])
                fint[jdofs] += fe

        for trFace in self.trNodes:
            self.trdofs = mesh.getDofIndices(trFace, self.T_doftypes)
            print("        fint = ", fint[self.trdofs])

        return max_hbw

    #-----------------------------------------------------------------------
    #   get_Int_Vector
    #-----------------------------------------------------------------------

    def get_Int_Vector(self, fint, disp, mesh):
        self.get_Matrix_0(None, fint, disp, mesh)

    #-----------------------------------------------------------------------
    #   get_Ext_Vector
    #-----------------------------------------------------------------------

    def get_Ext_Vector(self, fext, mesh):

        # Variables related to element on T mesh
        H = np.zeros((self.rank, self.ndof))
        eps = np.zeros((self.rank, self.rank))
        eps[0, 0] = self.strain[0]
        eps[1, 1] = self.strain[1]
        eps[0, 1] = eps[1, 0] = self.strain[2]/2

        # Loop over faces of trNodes
        for ix, trFace in enumerate(self.trNodes):

            u_corner = eps[ix, :]*self.dx_[ix]

            # Loop over indices of trFace
            for inod in range(len(trFace)-1):

                # Assemble H matrix
                connect = trFace[inod:inod+2]
                jdofs = mesh.getDofIndices(connect, self.T_doftypes)
                coords = mesh.getCoords(connect)
                w = self.bshape.getIntegrationWeights(coords)
                h = self.bshape.getShapeFunctions()

                # Vector to be assembled: fext[jdofs]
                fe = np.zeros(self.ndof)

                for ip in range(self.nIP):

                    # Assemble H matrix
                    H[0, 0] = H[1, 1] = h[ip][0]
                    H[0, 2] = H[1, 3] = h[ip][1]

                    # Assemble fe
                    fe = w[ip] * (H.transpose() @ u_corner)
                    fext[jdofs] += fe

        for trFace in self.trNodes:
            self.trdofs = mesh.getDofIndices(trFace, self.T_doftypes)
            print("        fext = ", fext[self.trdofs])

    #-----------------------------------------------------------------------
    #   get_Constraints
    #-----------------------------------------------------------------------

    def get_Constraints(self, cons, mesh):

        # Fix corner
        idofs = mesh.getDofIndices(self.corner0, self.U_doftypes)
        cons.addConstraints(idofs)
        
        # Voigt to tensor
        eps = np.zeros((self.rank, self.rank))
        eps[0, 0] = self.strain[0]
        eps[1, 1] = self.strain[1]
        eps[0, 1] = eps[1, 0] = self.strain[2]/2

        # Apply strain
        for ix in range(self.rank):
            for jx in range(self.rank):
                idof = mesh.getDofIndex(self.corner[ix], self.U_doftypes[jx])
                cons.addConstraint(idof, eps[ix, jx])

    #-----------------------------------------------------------------------
    #   takeAction
    #-----------------------------------------------------------------------

    def takeAction(self, action, mesh):
        if action == "plot_boundary":
            self.plotBoundary(mesh)

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

        for ix, trFace in enumerate(self.trNodes):
            coords = mesh.getCoords(trFace)
            coords[:, ix] -= self.dx_[ix]/10
            ax.plot(coords[:, 0], coords[:, 1],
                    marker='o', linewidth=0.3, markersize=3, color='blue')

        plt.show()

#===========================================================================
#   Example
#===========================================================================


if __name__ == '__main__':

    from properties import Properties
    from mesh import Mesh
    from models import ModelFactory
    from algebra import MatrixBuilder
    from constraints import Constraints
    np.set_printoptions(precision=4)

    # Props
    file = "Examples/rve.pro"
    props = Properties()
    props.parseFile(file)

    # Create mesh
    mesh = Mesh()
    mesh.initialize(props.getProps("input"))

    # Create model
    model = ModelFactory("model", props, mesh)

    # Create constraints
    cons = Constraints()
    cons.initialize(props, mesh)

    # Create matrix builder
    ndof = mesh.dofCount()
    mbuild = MatrixBuilder(ndof)

    # Create vectors
    fint = np.zeros(ndof)
    fext = np.zeros(ndof)
    disp = np.zeros(ndof)

    model.get_Constraints(cons, mesh)
    model.get_Ext_Vector(fext, mesh)
    model.get_Matrix_0(mbuild, fint, disp, mesh)
