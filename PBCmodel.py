# Import Standard Libraries
import scipy as np
import matplotlib.pyplot as plt
from copy import deepcopy

# Import Local Libraries
from models import Model
from shapes import Shape
from algebra import norm
from Voigt import voigt2TensorStrain

#===========================================================================
#   PBCmodel
#===========================================================================


class PBCmodel(Model):
    """ Periodic Boundary Conditions Model
        
    Instance Members:
        name = model name
        rank = number of dimensions
        type = model type ("Periodic")
        strain = prescribed strain rate
        factor = prescribed coarsening factor

        U_doftypes = displacement dof types
        T_doftypes = traction dof types
        bshape = boundary element shape
        nIP = number of integration points of bshape
        nnod = number of nodes associated with bshape
        ndof = number of degrees of freedom of bshape
        localrank = local rank (ndim) of bshape

        bndNodes = boundary nodes = [xmin, xmax, ymin, ymax]
        trNodes = traction nodes = [xmin, ymin]
        corner0 = corner at xmin & ymin
        corner = [cornerx, cornery]

        dx_0 = smallest element dimensions
        box = specimen coordinates
        dx = specimen dimensions

    Public Methods:
        PBCModel(name, conf, props, globdat)
        takeAction(action, globdat)

    Private Methods:
        __boundingBox(mesh)
        __setTolerances()
        __findBndNodes(mesh)
        __sortBndNodes(mesh)
        __sortBndFace(mesh, bndFace, index)
        __findCornerNodes()
        __findSmallestElement(mesh)
        __createTractionMesh(mesh)
        __coarsenMesh(mesh, trFace, index)
        __get_Matrix_0(mbuild, fint, disp, mesh)
        __getTractionMeshNodes(mesh, x, face)
        __get_Constraints(cons, mesh)
        __plotMeshAndBoundary(mesh)
        __plotBoundary(mesh)
        __advance(i)
    """

    # Public:

    #-----------------------------------------------------------------------
    #   constructor
    #-----------------------------------------------------------------------

    def __init__(self, name, conf, props, globdat):
        self.name = name
        myConf = conf.makeProps(name)
        myProps = props.getProps(name)

        self.type = myProps.get("type", "Periodic")
        self.strain = myProps.find("strainRate", None)
        self.factor = myProps.get("coarsenFactor",1.0)

        myConf.set("type", self.type)
        myConf.set("strainRate", self.strain)
        myConf.set("coarsenFactor", self.factor)

        mesh = globdat.get("mesh")
        self.rank = mesh.rank

    #-----------------------------------------------------------------------
    #   initialize
    #-----------------------------------------------------------------------

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
        self.bshape = Shape.shapeFactory(myConf, myProps)
        self.localrank = self.bshape.ndim
        self.nnod = self.bshape.nnod
        self.nIP = self.bshape.nIP
        self.ndof = self.nnod*self.rank
        if self.localrank != self.rank-1:
            msg = "Shape ndim = {}. Should be {}".format(
                self.localrank, self.rank-1)
            raise ValueError(msg)
            
    #-----------------------------------------------------------------------
    #   takeAction
    #-----------------------------------------------------------------------

    def takeAction(self, action, globdat):
        if action == "GET_MATRIX_0" or action == "GET_INT_VECTOR":
            mesh = globdat.get("mesh")
            mbuild = globdat.get("mbuild")
            fint = globdat.get("fint")
            disp = globdat.get("solu")
            self.__get_Matrix_0(mbuild, fint, disp, mesh)
            return True
        elif action == "GET_CONSTRAINTS":
            mesh = globdat.get("mesh")
            cons = globdat.get("cons")
            self.__get_Constraints(cons, mesh)
            return True
        elif action == "PLOT_BOUNDARY":
            mesh = globdat.get("mesh")
            self.__plotBoundary(mesh)
            return True
        elif action == "PLOT_MESH":
            mesh = globdat.get("mesh")
            self.__plotMeshAndBoundary(mesh)
            return True
        elif action == "ADVANCE":
            self.__advance(globdat.i)
            return True
        else:
            return False

    # Private:

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
        self.dx = np.zeros(self.rank)

        # Find specimen coordinates
        for coord in mesh.coords:
            for ix in range(self.rank):
                if coord[ix] < self.box_[2*ix]:
                    self.box_[2*ix] = coord[ix]
                if coord[ix] > self.box_[2*ix+1]:
                    self.box_[2*ix+1] = coord[ix]

        # Find specimen dimensions
        for ix in range(self.rank):
            self.dx[ix] = self.box_[2*ix + 1] - self.box_[2*ix]

        # Print to verify
        print("        box = ", self.box_)
        print("        dx = ", self.dx)

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
    #   __findCornerNodes
    #-----------------------------------------------------------------------

    def __findCornerNodes(self):
        """ Finds the intersection of each bndFace """

        self.corner = []
        if self.rank == 2:
            self.corner0 = self.bndNodes[0][0]
            self.corner.append(self.bndNodes[1][0])
            self.corner.append(self.bndNodes[3][0])
            print("        corner0 = ", self.corner0)
            print("        cornerx = ", self.corner[0])
            print("        cornery = ", self.corner[1])
        elif self.rank == 3:
            raise NotImplementedError(" Not yet implemented. ")
            
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
    #   __findSmallestElement
    #-----------------------------------------------------------------------

    def __findSmallestElement(self, mesh):
        """ Finds the smallest element dimension on each bndFace """

        # smallest element dimensions
        self.dx_0 = deepcopy(self.dx)

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
                if dx < self.dx_0[index]:
                    self.dx_0[index] = dx

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
                coord[ix] = self.box_[2*ix + 1]
                trFace.append(mesh.addNode(coord))

            # Loop over indices of jnodes
            for jnod in jnodes:
                coord = mesh.getCoords(jnod)
                trFace.append(mesh.addNode(coord))

            # Get correct index
            index = 1 if (ix == 0) else 0

            # Sorting trFace (works only for 2D)
            self.__sortBndFace(mesh, trFace, index)

            # Coarsen trFace (works only for 2D)
            self.__coarsenMesh(mesh, trFace, index)

            # Add dofs to traction mesh
            mesh.addDofs(trFace, self.T_doftypes)
            
            # Print to verify
            print("        trNodes[{}] = {}".format(ix, self.trNodes[ix]))


    #-----------------------------------------------------------------------
    #   __coarsenMesh
    #-----------------------------------------------------------------------

    def __coarsenMesh(self, mesh, trFace, index):
        """ Coarsens the traction mesh """

        cn = mesh.getCoords(trFace[-1])
        dx = (self.dx_0[0]+self.dx_0[1])/(2*self.factor)

        # Loop over indices of trFace:
        for inod in range(len(trFace)-1):

            # Get nodal coordinates
            c0 = mesh.getCoords(trFace[inod])
            c1 = mesh.getCoords(trFace[inod+1])

            # Delete indices until c1 - c0 > dx
            while norm(c1 - c0) < min(dx,self.dx[index]):
                # Delete current index
                del trFace[inod+1]
                # Assign next node coords to c1
                c1 = mesh.getCoords(trFace[inod+1])

            # Check distance to last node
            if norm(cn - c1) < dx:
                # Delete all indices up to but not including the last one
                del trFace[inod+1:-1]
                return

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
            
            raise RuntimeError(" No connect found. ")

        elif self.rank == 3:
            raise NotImplementedError(" Not yet implemented. ")

    #-----------------------------------------------------------------------
    #   __get_Matrix_0
    #-----------------------------------------------------------------------

    def __get_Matrix_0(self, mbuild, fint, disp, mesh):
        """ Augments mbuild and fint with Ke, KeT, H, etc. """

        # Loop over faces of bndNodes
        for face, bndFace in enumerate(self.bndNodes):

            # Loop over indices of bndFace
            for inod in range(len(bndFace)-1):

                # Get idofs and w from U-mesh
                connect = bndFace[inod:inod+2]
                coords = mesh.getCoords(connect)
                X = self.bshape.getGlobalPoints(coords)
                w = self.bshape.getIntegrationWeights(coords)
                idofs = mesh.getDofIndices(connect, self.U_doftypes)

                for ip in range(self.nIP):

                    # Get jdofs from T-mesh
                    connect = self.__getTractionMeshNodes(mesh, X[ip], face)
                    jdofs = mesh.getDofIndices(connect, self.T_doftypes)
                    coords = mesh.getCoords(connect)

                    # Get N-matrix from U-mesh
                    N = self.bshape.getNmatrix(ip, ndim=self.rank)

                    # Get H-matrix from T-mesh
                    xi = self.bshape.getLocalPoint(X[ip], coords)
                    H = self.bshape.evalNmatrix(xi, ndim=self.rank)

                    # Assemble Ke
                    if face == 0 or face == 2 or face == 4:
                        # Negative face
                        Ke = + w[ip] * N.transpose() @ H
                    elif face == 1 or face == 3 or face == 5:
                        # Possitive face
                        Ke = - w[ip] * N.transpose() @ H

                    # Add Ke and KeT to mbuild
                    KeT = Ke.transpose()
                    mbuild.addBlock(idofs, jdofs, Ke)
                    mbuild.addBlock(jdofs, idofs, KeT)

                    # Assemble U-mesh fe
                    fe = Ke.dot(disp[jdofs])
                    fint[idofs] += fe

                    # Assemble T-mesh fe
                    fe = KeT.dot(disp[idofs])
                    fint[jdofs] += fe

        # Variables related to corner displacements
        kdofs = mesh.getDofIndices(self.corner0, self.U_doftypes)

        # Loop over faces of trNodes
        for ix, trFace in enumerate(self.trNodes):

            idofs = mesh.getDofIndices(self.corner[ix], self.U_doftypes)
            u_corner = disp[idofs]

            # Loop over indices of trFace
            for inod in range(len(trFace)-1):

                # Assemble H matrix
                connect = trFace[inod:inod+2]
                jdofs = mesh.getDofIndices(connect, self.T_doftypes)
                coords = mesh.getCoords(connect)
                w = self.bshape.getIntegrationWeights(coords)

                for ip in range(self.nIP):

                    # Assemble H matrix
                    H = w[ip] * self.bshape.getNmatrix(ip, ndim=self.rank)
                    Ht = H.transpose()

                    # Add H matrix to mbuild
                    mbuild.addBlock(idofs, jdofs, H)
                    mbuild.addBlock(jdofs, idofs, Ht)

                    # Assemble U-mesh fe
                    fint[idofs] += H @ disp[jdofs]
                    fint[kdofs] -= H @ disp[jdofs]

                    # Assemble T-mesh fe
                    fint[jdofs] += Ht @ u_corner

    #-----------------------------------------------------------------------
    #   __get_Constraints
    #-----------------------------------------------------------------------

    def __get_Constraints(self, cons, mesh):

        # Fix corner
        idofs = mesh.getDofIndices(self.corner0, self.U_doftypes)
        cons.addConstraints(idofs)
        
        # Voigt to tensor
        if self.strain is not None:
            eps = voigt2TensorStrain(self.strain)

            # Apply strain
            for ix in range(self.rank):
                for jx in range(self.rank):
                    idof = mesh.getDofIndex(self.corner[ix], self.U_doftypes[jx])
                    cons.addConstraint(idof, eps[ix, jx]*self.dx[ix])

    #-----------------------------------------------------------------------
    #   __plotBoundary
    #-----------------------------------------------------------------------

    def __plotBoundary(self, mesh):
        """ Plots the boundary nodes """

        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111)

        for bndFace in self.bndNodes:
            coords = mesh.getCoords(bndFace)
            ax.plot(coords[:, 0], coords[:, 1],
                    marker='o', linewidth=0.3, markersize=3, color='k')

        for ix, trFace in enumerate(self.trNodes):
            coords = mesh.getCoords(trFace)
            coords[:, ix] += self.dx[ix]/10
            ax.plot(coords[:, 0], coords[:, 1],
                    marker='o', linewidth=0.3, markersize=3, color='blue')

        plt.show()

    #-----------------------------------------------------------------------
    #   __plotMesh
    #-----------------------------------------------------------------------

    def __plotMeshAndBoundary(self, mesh):
        """ Plots the mesh and boundary nodes """

        ax = mesh.plotMesh()
        for bndFace in self.bndNodes:
            coords = mesh.getCoords(bndFace)
            ax.plot(coords[:, 0], coords[:, 1],
                    marker='o', linewidth=0.3, markersize=3, color='k')

        for ix, trFace in enumerate(self.trNodes):
            coords = mesh.getCoords(trFace)
            coords[:, ix] += self.dx[ix]/10
            ax.plot(coords[:, 0], coords[:, 1],
                    marker='o', linewidth=0.3, markersize=3, color='blue')

        plt.show()

    #-----------------------------------------------------------------------
    #   __advance
    #-----------------------------------------------------------------------

    def __advance(self, i):
        if self.strain is not None:
            pass
