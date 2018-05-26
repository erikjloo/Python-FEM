# Import Standard Libraries
import re
import scipy as np

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
        SolidModel(name, props, mesh)
        get_Matrix_0(mbuild, fint, disp, mesh)
        get_Ext_Vector(fext, mesh)
        get_Int_Vector(fint, disp, mesh)
        get_Constraints(cons, mesh)
        takeAction(action, mesh)

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
        if self.rank < 1 or self.rank > 3:
            msg = "Rank = {}. Should be 1, 2 or 3".format(self.rank)
            raise ValueError(msg)

        # Get element group
        if mesh.doElemGroups is True:
            gmsh_group = props.get("elements")
            key = int(re.search(r'\d+', gmsh_group).group())
            group_name = mesh.groupNames[key]
            print("    Obtaining elements from {}".format(group_name))
            idx = mesh.groupNames.keys().index(key)
            self.ielements = mesh.groups[idx]
            print("    Elements in mesh.groups[{}]".format(idx))
        else:
            group_name = next(iter(mesh.groupNames.values()))
            self.ielements = mesh.groups[0]
            print("    Obtaining elements from {}".format(group_name))

        # Add types
        types = ['u', 'v', 'w']
        self.types = [ types[x] for x in range(self.rank)]
        mesh.addTypes(self.types)

        # Add dofs
        self.inodes = mesh.getNodeIndices(self.ielements)
        mesh.addDofs(self.inodes, self.types)

        # Add thickness (2D)
        if self.rank == 2:
            self.t = props.get("thickness")

        # Create element
        self.shape = ShapeFactory(props)
        localrank = self.shape.ndim
        if localrank != self.rank:
            msg = "Shape ndim = {}. Should be {}".format(localrank,self.rank)
            raise ValueError(msg)

        # Create material
        self.mat = MaterialFactory(props)
        
        
    #-----------------------------------------------------------------------
    #   get_Matrix_0
    #-----------------------------------------------------------------------

    def get_Matrix_0(self, mbuild, fint, disp, mesh):
        """ Input & Output: mbuild = MatrixBuilder = Ksys
                            fint = internal force vector """

        max_hbw = 0
        print(disp)
        # Iterate over elements assigned to model
        for iele in self.ielements:

            # Get element nodes, coordinates, dofs and displacements
            inodes = mesh.getNodes(iele)
            coords = mesh.getCoords(inodes)
            idofs = mesh.getDofIndices(inodes, self.types)
            hbw = max(idofs) - min(idofs)
            max_hbw = hbw if hbw > max_hbw else max_hbw
            ele_disp = disp[idofs]

            # Initialize element stiffness matrix and internal force vector
            ndof = len(idofs)
            kele = np.zeros((ndof, ndof))
            fe = np.zeros(ndof)

            for ip in range(self.shape.nIP):

                # Get strain, B matrix and weight
                [strain, B, w] = self.shape.getStrain(coords, ele_disp, ip)
                w = w*self.t if self.rank == 2 else w
                print(iele)
                print(ele_disp)
                # print(strain)
                # Get tangent stiffness matrix D
                [stress, D] = self.mat.getStress(strain)

                # Compute element stiffness matrix (kele)
                kele += w * (B.transpose() @ D @ B)

                # Compute the element internal force vector
                fe += w * (B.transpose() @ stress)

            # Add kele to the global stiffness matrix (Ksys):
            mbuild.addBlock(idofs, idofs, kele)

            # Add fint to the global force vector (Fint):
            fint[idofs] += fe
            
        return max_hbw

    #-----------------------------------------------------------------------
    #   get_Int_Vector
    #-----------------------------------------------------------------------

    def get_Int_Vector(self, fint, disp, mesh):
        """ Input & Output: fint = internal force vector """

        # Iterate over elements assigned to model
        for iele in self.ielements:

            # Get element nodes, coordinates, dofs and displacements
            inodes = mesh.getNodes(iele)
            coords = mesh.getCoords(inodes)
            idofs = mesh.getDofIndices(inodes, self.types)
            ele_disp = disp[idofs]
            
            # Initialize element internal force vector
            fe = np.zeros(len(idofs))

            for ip in range(self.shape.nIP):

                # Get strain, B matrix and weight
                [strain, B, w] = self.shape.getStrain(coords, ele_disp, ip)
                w = w*self.t if self.rank == 2 else w

                # Get tangent stiffness matrix D
                [stress, _] = self.mat.getStress(strain)
                
                # Compute the element internal force vector
                fe += w * (B.transpose() @ stress)

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

    def get_Constraints(self, cons, mesh):
        pass

    #-----------------------------------------------------------------------
    #   takeAction
    #-----------------------------------------------------------------------

    def takeAction(self, action, *args):
        pass
