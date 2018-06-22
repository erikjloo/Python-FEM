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
        name = model name
        type = model type ("Solid")
        rank = number of dimensions
        group = element group
        
        types = displacement dof types
        ielements = element indices
        inodes = node indices 
        
        shape = element shape
        mat = material

    Public Methods:
        SolidModel(name, props, mesh)
        takeAction(action, globdat)
    
    Private Methods:
        __get_Matrix_0(mbuild, fint, disp, mesh)
        __get_Int_Vector(fint, disp, mesh)

    """

    # Public:

    #-----------------------------------------------------------------------
    #   constructor
    #-----------------------------------------------------------------------

    def __init__(self, name, props, conf, mesh):

        self.name = name
        self.rank = mesh.rank
        myConf = conf.makeProps(name)
        myProps = props.getProps(name)

        self.type = myProps.get("type","Solid")
        self.group = myProps.get("elements","All")

        myConf.set("type", self.type)
        myConf.set("elements", self.group)

        if self.group != "All":
            key = int(re.search(r'\d+', self.group).group())
            group_name = mesh.groupNames[key]
            print("    Obtaining elements from {}".format(group_name))
            idx = mesh.groupNames.keys().index(key)
            self.ielements = mesh.groups[idx]
            print("    Elements in mesh.groups[{}]".format(idx))
        else:
            group_name = next(iter(mesh.groupNames.values()))
            self.ielements = mesh.groups[0]
            print("    Obtaining elements from {}".format(group_name))
        
    #-----------------------------------------------------------------------
    #   initialize
    #-----------------------------------------------------------------------

        # Add types
        types = ['u', 'v', 'w']
        self.types = [ types[x] for x in range(self.rank)]
        mesh.addTypes(self.types)

        # Add dofs
        self.inodes = mesh.getNodeIndices(self.ielements)
        mesh.addDofs(self.inodes, self.types)

        # Add thickness (2D)
        if self.rank == 2:
            self.t = myProps.get("thickness",1.0)
            myConf.set("thickness",self.t)

        # Create element
        self.shape = ShapeFactory(myProps, myConf)
        localrank = self.shape.ndim
        if localrank != self.rank:
            msg = "Shape ndim = {}. Should be {}".format(localrank,self.rank)
            raise ValueError(msg)

        # Create material
        self.mat = MaterialFactory(myProps, myConf)
        
    #-----------------------------------------------------------------------
    #   takeAction
    #-----------------------------------------------------------------------

    def takeAction(self, action, globdat):
        if action == "GET_MATRIX_0":
            self.__get_Matrix_0(globdat.mbuild, globdat.fint, globdat.disp, globdat.mesh)
            return True
        elif action == "GET_INT_VECTOR":
            self.__get_Int_Vector(globdat.fint, globdat.disp, globdat.mesh)
            return True
        else:
            return False

    #-----------------------------------------------------------------------
    #   get_Matrix_0
    #-----------------------------------------------------------------------

    def __get_Matrix_0(self, mbuild, fint, disp, mesh):
        """ Input & Output: mbuild = MatrixBuilder = Ksys
                            fint = internal force vector """

        max_hbw = 0
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

                # Get tangent stiffness matrix D
                [stress, D] = self.mat.getStress(strain)

                # Compute element stiffness matrix (kele)
                kele += w * (B.transpose() @ D @ B)

                # Compute the element internal force vector
                fe += w * (B.transpose() @ stress)

            # Add kele to the global stiffness matrix (Ksys):
            mbuild.addBlock(idofs, idofs, kele)

            # Add fint to the global force vector (Fint):
            fint[idofs] += kele @ disp[idofs]
            
        mbuild.hbw = max_hbw if mbuild.hbw < max_hbw else mbuild.hbw

    #-----------------------------------------------------------------------
    #   get_Int_Vector
    #-----------------------------------------------------------------------

    def __get_Int_Vector(self, fint, disp, mesh):
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

