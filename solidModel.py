# Import Standard Libraries
import scipy as np

#Import Local Libraries
from mesh import Mesh
from dofspace import DofSpace
from properties import Properties
from elements import ElementType
from algebra import MatrixBuilder
from shapes import Line2, Tri3, Quad4


#===========================================================================
#   SolidModel
#===========================================================================


class SolidModel(Mesh, DofSpace):
    """ Solid Model

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

        props = list of element type and physical group of each element
        Phys = names of physical groups
        nPhys = number of physical groups

        dofspace = array of dof indices (idofs)
        types = list of dof type names
        idof = last dof index = ndof - 1
        ndof = number of dofs

        rank = number of dimensions
        
        shape = element shape

        material = ?

    Public Methods:
        SolidModel()

    Private Methods:
        __verifyShape()
        __getBmatrix()
        __getStrain()
        __getStress()
    """

    # Static:
    __rank_error__ = "Rank has to be 1, 2 or 3!"

    # Public:

    #-----------------------------------------------------------------------
    #   Constructor
    #-----------------------------------------------------------------------

    def __init__(self, name, props, rank=2):

        # myProps = props.getProperties(name)
        # myConf = props.makeProperties(name)

        # Call the Mesh constructor
        path = props.get("userInput.mesh.file")
        Mesh.__init__(self)
        self.readMesh(path)
        self.rank = rank

        # Call the DofSpace constructor
        DofSpace.__init__(self, self.nnod, self.rank)


    def initialize(self):

        # Add types and shape
        if self.rank == 1:
            self.addType("u")
            self.shape = Line2()
        elif self.rank == 2:
            self.addTypes(["u", "v"])
            self.shape = Tri3()
        elif self.rank == 3:
            self.addTypes(["u", "v", "w"])
            self.shape = Quad4()
        else:
            raise ValueError(self.__rank_error__)

        # Add dofs
        self.addDofs(range(self.nnod), self.types)

    def get_Matrix_0(self, mbuild, F_int):
        """ Input & Output: mbuild = MatrixBuilder = Ksys
                            F_int = internal force vector """

        # Iterate over elements assigned to model
        for iele in range(self.nele):

            # Get element nodes
            inodes = self.getNodes(iele)

            # Get nodal coordintates
            coords = self.getCoords(inodes)

            # Get the element degrees of freedom
            idofs = self.getDofIndices(inodes, self.types)

            # pylint: disable = unbalanced-tuple-unpacking

            # get nodal displacements

            # multiply w[ip] times thickness

            # element stiffness matrix
            ndof = len(idofs)
            kele = np.zeros((ndof, ndof))

            E, v = 29000, 0.3
            la = v*E/((1+v)*(1-2*v))
            mu = E/(2*(1+v))
            D = np.array([[la+2*mu, la, 0], [la, la+2*mu, 0], [0, 0, mu]])

            # self.__verifyShape(etype)
            B, w = self.__getBmatrix(coords)

            for ip in range(self.shape.nIP):

                # strain, eledisp = shape.getstrain( ip, ie)

                # Get tangent stiffness matrix:
                # material.update(stress, D, strain, ip)

                # Compute element stiffness matrix (kele):
                kele += (B[ip].transpose() @ D @ B[ip])*w[ip]
            
                # Compute the element force vector (fint):
                # fint += w[ip] * B[ip].transpose @ stress

            # Add k to the global stiffness matrix (Ksys):
            mbuild.addBlock(idofs, idofs, kele)

            # Add fint to the global force vector (Fint):
            # Fint[idofs] += fint

        return mbuild, F_int
        # output:
        # u, strain, sigma, f_int, K

    def __verifyShape(self, etype):
        """ Verifies if a new shape type is necessary """
        if etype == ElementType.line2 and not isinstance(self.shape, Line2):
            self.shape = Line2()
        elif etype == ElementType.tri3 and not isinstance(self.shape, Tri3):
            self.shape = Tri3()
        elif etype == ElementType.quad4 and not isinstance(self.shape, Quad4):
            self.shape = Quad4()

    def __getBmatrix(self, coords):
        """ Returns the B and w given the element nodal coordinates """
        if self.rank == 1:
            B, w = self.shape.getGlobalGradients(coords[:, 0])
        elif self.rank == 2:
            B, w = self.shape.getBmatrix(coords[:, 0:2])
        elif self.rank == 3:
            B, w = self.shape.getBmatrix(coords)
        return B, w


#===========================================================================
#   Example
#===========================================================================


if __name__ == '__main__':

    file = "Examples/square.pro"
    props = Properties()
    props.parseFile(file)

    model = SolidModel("model.matrix", props, rank=2)
    model.initialize()

    ndof = model.dofCount()
    f_int = np.zeros(ndof)
    mbuild = MatrixBuilder(ndof)
    mbuild, f_int = model.get_Matrix_0(mbuild, f_int)
    mbuild.print()
