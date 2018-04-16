# Import Standard Libraries
import scipy as np

#Import Local Libraries
from mesh import Mesh
from dofspace import DofSpace
from properties import Properties, ElementType
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
        __rank_error__ = "Rank has to be between 1 and 3"
        
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

        shape = element shape (may change if more than one element type)
        rank = number of dimensions

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
    __rank_error__ = "Rank has to be between 1 and 3"

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

        self.readMesh(path)
        self.rank = rank
        self.dofspace = np.empty((self.nnod, rank))
        self.dofspace[:] = np.nan
        # Add types and shape
        if rank == 1:
            self.addType("u")
            self.shape = Line2()
        elif rank == 2:
            self.addTypes(["u", "v"])
            self.shape = Tri3()
        elif rank == 3:
            self.addTypes(["u", "v", "w"])
            self.shape = Quad4()
        else:
            raise ValueError(self.__rank_error__)

        # Add dofs
        self.addDofs(range(self.nnod), self.types)

        self.props = Properties(self.nele, self.props)
        self.props.addMaterial("Matrix", E=9000, v=0.3)
        self.props.addMaterial("Fibers", E=29000, v=0.2)

    def assemble(self, mbuild, F_int):
        """ Input & Output: mbuild = MatrixBuilder = Ksys
                            F_int = internal force vector """

        # Iterate over elements assigned to model
        for iele, inodes in enumerate(self.connectivity):

            # Get nodal coordintates
            coords = self.getCoords(inodes)

            # Get the element degrees of freedom
            idofs = self.getDofIndices(inodes, self.types)

            # pylint: disable = unbalanced-tuple-unpacking
            [etype, mat] = self.props.getProperties(iele, ["etype", "mat"])

            # get nodal displacements

            # multiply w[ip] times thickness

            # element stiffness matrix
            ndof = len(idofs)
            kele = np.zeros((ndof, ndof))

            E, v = mat.E, mat.v
            la = v*E/((1+v)*(1-2*v))
            mu = E/(2*(1+v))
            D = np.array([[la+2*mu, la, 0], [la, la+2*mu, 0], [0, 0, mu]])

            self.__verifyShape(etype)
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

    model = SolidModel("Examples/rve.msh", rank=2)

    ndof = model.dofCount()
    F_int = np.zeros(ndof)
    mbuild = MatrixBuilder(ndof)
    mbuild = model.assemble(mbuild, F_int)
