# Import Standard Libraries
import scipy as np
import re

# Import Local Libraries
from models import Model
from shapes import ShapeFactory
from materials import MaterialFactory


#===========================================================================
#   SolidModel
#===========================================================================


class SolidModel(Model):
    """ Solid Model
        
    Instance Members:
        ielements = element indices
        rank = number of dimensions
        types = displacement dof types

        shape = element shape
        nIP = number of integration points of shape
        nnod = number of nodes of shape
        localrank = local rank of shape

    Public Methods:
        SolidModel()
        get_Matrix_0(mesh, mbuild, fint)
        get_Ext_Vector(fext)
        get_Constraints(mesh, constraints)

    Private Methods:
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

        # Add types
        types = ['u', 'v', 'w']
        self.types = [ types[x] for x in range(self.rank)]
        mesh.addTypes(self.types)

        # Add dofs
        inodes = mesh.getNodeIndices(self.ielements)
        mesh.addDofs(inodes, self.types)

        # Create element
        self.shape = ShapeFactory(props)
        self.nIP = self.shape.nIP
        self.nnod = self.shape.nnod
        self.localrank = self.shape.ndim

        # Create material
        self.mat = MaterialFactory(props)
        self.nIP = self.shape.nIP
        self.nnod = self.shape.nnod
        self.localrank = self.shape.ndim

    #-----------------------------------------------------------------------
    #   get_Matrix_0
    #-----------------------------------------------------------------------

    def get_Matrix_0(self, mbuild, fint, disp, mesh):
        """ Input & Output: mbuild = MatrixBuilder = Ksys
                            fint = internal force vector """

        # Iterate over elements assigned to model
        for iele in self.ielements:

            # Get element nodes, coordinates, dofs and displacements
            inodes = mesh.getNodes(iele)
            coords = mesh.getCoords(inodes)
            idofs = mesh.getDofIndices(inodes, self.types)
            ele_disp = disp[idofs]

            # Initialize element stiffness matrix and internal force vector
            ndof = len(idofs)
            kele = np.zeros((ndof, ndof))
            fe = np.zeros(ndof)

            for ip in range(self.shape.nIP):

                # Get strain, B matrix and weight
                [strain, B, w] = self.shape.getStrain(coords, ele_disp, ip)

                # Get tangent stiffness matrix D
                [stress, D] = self.mat.getStress(strain[ip])

                # Compute element stiffness matrix (kele)
                kele += w * (B.transpose() @ D @ B)
            
                # Compute the element internal force vector
                fe += w * (B.transpose @ stress)

            # Add kele to the global stiffness matrix (Ksys):
            mbuild.addBlock(idofs, idofs, kele)

            # Add fint to the global force vector (Fint):
            fint[idofs] += fe

    #-----------------------------------------------------------------------
    #   get_Ext_Vector
    #-----------------------------------------------------------------------

    def get_Ext_Vector(self, fext, mesh):
        pass

    #-----------------------------------------------------------------------
    #   get_Constraints
    #-----------------------------------------------------------------------

    def get_Constraints(self, mesh, constraints):
        pass

    #-----------------------------------------------------------------------
    #   takeAction
    #-----------------------------------------------------------------------

    def takeAction(self, action, *args):
        pass

#===========================================================================
#   Example
#===========================================================================


if __name__ == '__main__':

    from properties import Properties
    from algebra import MatrixBuilder
    from models import ModelFactory
    from nonlin import multistep
    from mesh import Mesh

    # Initialization
    file = "Examples/semicircle.pro"
    props = Properties()
    props.parseFile(file)

    # Mesh
    mesh = Mesh()
    mesh.initialize(props.getProps("input.mesh"))

    # Model
    model = ModelFactory("model", props, mesh)

    ndof = mesh.dofCount()
    mbuild = MatrixBuilder(ndof)
    fint = np.zeros(ndof)
    disp = np.zeros(ndof)

    model.get_Matrix_0(mbuild, fint, disp, mesh) 
