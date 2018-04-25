# Import Standard Libraries
import scipy as np
import re

#Import Local Libraries
from models import Model
from algebra import MatrixBuilder
from shapes import ShapeFactory


#===========================================================================
#   SolidModel
#===========================================================================


class SolidModel(Model):
    """ Solid Model
        
    Instance Members:
        ielements = element indices
        rank = number of dimensions

        shape = element shape
        nIP = number of integration points of shape
        nnod = number of nodes of shape
        localrank = local rank of shape

        types = displacement dof types

        
    Public Methods:
        SolidModel()

    Private Methods:
        __verifyShape()
        __getBmatrix()
        __getStrain()
        __getStress()
    """

    # Public:

    #-----------------------------------------------------------------------
    #   initialize
    #-----------------------------------------------------------------------

    def __init__(self, name, props, mesh):

        # Model name
        self.name = name
        self.rank = mesh.rank

        # Get element group
        gmsh_group = props.get("elements")
        group = int(re.search(r'\d+', gmsh_group).group())
        self.ielements = mesh.groups[group]

        # Create element
        self.shape = ShapeFactory(props)
        self.nIP = self.shape.nIP
        self.nnod = self.shape.nnod
        self.localrank = self.shape.ndim

        # Add types
        types = ['u', 'v', 'w']
        self.types = [ types[x] for x in range(self.rank)]
        mesh.addTypes(self.types)

        # Add dofs
        inodes = mesh.getNodeIndices(self.ielements)
        mesh.addDofs(inodes, self.types)

    #-----------------------------------------------------------------------
    #   get_Matrix_0
    #-----------------------------------------------------------------------

    def get_Matrix_0(self, mesh, mbuild, f_int):
        """ Input & Output: mbuild = MatrixBuilder = Ksys
                            F_int = internal force vector """

        # Iterate over elements assigned to model
        for iele in self.ielements:

            # Get element nodes
            inodes = mesh.getNodes(iele)

            # Get nodal coordintates
            coords = mesh.getCoords(inodes)

            # Get the element degrees of freedom
            idofs = mesh.getDofIndices(inodes, self.types)

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

        # output:
        # u, strain, sigma, f_int, K

    #-----------------------------------------------------------------------
    #   get_Ext_Vector
    #-----------------------------------------------------------------------

    def get_Ext_Vector(self, f_ext):
        pass

    #-----------------------------------------------------------------------
    #   get_Constraints
    #-----------------------------------------------------------------------

    def get_Constraints(self, mesh, constraints):
        pass

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

    from properties import Properties
    from models import ModelFactory
    from mesh import Mesh

    file = "Examples/rve.pro"
    props = Properties()
    props.parseFile(file)

    mesh = Mesh()
    mesh.initialize(props, rank=2)

    model = ModelFactory("model", props, mesh)

    ndof = mesh.dofCount()
    f_int = np.zeros(ndof)
    mbuild = MatrixBuilder(ndof)

    model.get_Matrix_0(mesh, mbuild, f_int)
    mbuild.print()
