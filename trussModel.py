# Import Standard Libraries
import re
import scipy as np
from copy import deepcopy

# Import Local Libraries
from models import Model
from algebra import norm, gram_schmidt

#===========================================================================
#   TrussModel
#===========================================================================


class TrussModel(Model):
    """ Truss Model
        
    Instance Members:
        name = model name
        type = model type ("Truss")
        rank = number of dimensions
        group = element group
        
        types = displacement dof types
        ielements = element indices
        inodes = node indices
        
        E = Young's modulus
        A = cross-section area

    Public Methods:
        TrussModel(name, conf, props, globdat)
        takeAction(action, globdat)
    
    Private Methods:
        __get_Matrix_0(mbuild, fint, du, mesh)
        __get_Int_Vector(fint, du, mesh)

    """

    # Public:

    #-----------------------------------------------------------------------
    #   constructor
    #-----------------------------------------------------------------------

    def __init__(self, name, conf, props, globdat):
        self.name = name
        myConf = conf.makeProps(name)
        myProps = props.getProps(name)

        self.type = myProps.get("type", "Truss")
        self.group = myProps.get("elements", "All")

        myConf.set("type", self.type)
        myConf.set("elements", self.group)

        mesh = globdat.get("mesh")
        self.rank = mesh.rank

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

        nele = len(self.ielements)
        self.stress0 = np.zeros(nele)
        self.stress = np.zeros(nele)
        self.strain0 = np.zeros(nele)
        self.strain = np.zeros(nele)

        # Add types
        types = ['u', 'v', 'w']
        self.types = [types[x] for x in range(self.rank)]
        mesh.addTypes(self.types)

        # Add dofs
        self.inodes = mesh.getNodeIndices(self.ielements)
        mesh.addDofs(self.inodes, self.types)

        # Create element
        self.A = myProps.get("area", 1.0)
        self.E = myProps.get("young", 1.0)
        myConf.set("area", self.A)
        myConf.set("young", self.E)

        self.L = [self.__length(iele, mesh) for iele in self.ielements]

    #-----------------------------------------------------------------------
    #   takeAction
    #-----------------------------------------------------------------------

    def takeAction(self, action, globdat):
        if action == "GET_MATRIX_0":
            mbuild = globdat.get("mbuild")
            fint = globdat.get("fint")
            du = globdat.get("du")
            mesh = globdat.get("mesh")
            self.__get_Matrix_0(mbuild, fint, du, mesh)
            return True
        elif action == "GET_INT_VECTOR":
            fint = globdat.get("fint")
            du = globdat.get("du")
            mesh = globdat.get("mesh")
            self.__get_Matrix_0(None, fint, du, mesh)
            return True
        else:
            return False

    #-----------------------------------------------------------------------
    #   get_Matrix_0
    #-----------------------------------------------------------------------

    def __get_Matrix_0(self, mbuild, fint, du, mesh):
        EA = self.E*self.A
        mesh.updateGeometry(du)

        for iele in self.ielements:

            Lo = self.L[iele]
            inodes = mesh.getNodes(iele)
            coords = mesh.getCoords(inodes)
            idofs = mesh.getDofIndices(inodes, self.types)

            Gamma, L = self.__transformation(coords)

            K_l = np.zeros((2*self.rank, 2*self.rank))
            K_l[[0, self.rank], [0, self.rank]] = EA/Lo
            K_l[[0, self.rank], [self.rank, 0]] = -EA/Lo
            K_l = Gamma.transpose() @ K_l @ Gamma

            stress = self.E*(L-Lo)/Lo

            K_nl = np.eye(4)*self.A*stress/Lo
            K_nl[[0, 1, 2, 3], [2, 3, 0, 1]] = -self.A*stress/Lo

            if mbuild is not None:
                mbuild.addBlock(idofs, idofs, K_l + K_nl)
            fe = self.A*np.array([-stress,0,stress,0])
            fint[idofs] += Gamma.transpose().dot(fe)

    def __length(self, iele, mesh):
        inodes = mesh.getNodes(iele)
        coords = mesh.getCoords(inodes)
        return norm(coords[1, :] - coords[0, :])

    def __transformation(self, coords):
        Gamma = np.zeros((2*self.rank, 2*self.rank))
        if self.rank == 2:
            dx_dy = coords[1, :] - coords[0, :] # [x2-x1 y2-y1]
            i_bar = dx_dy/norm(dx_dy) 
            j_bar = i_bar.dot([[0, 1], [-1, 0]]) # Rotate i_bar by 90 degrees
            Gamma[0, 0:2] = Gamma[2, 2:4] = i_bar
            Gamma[1, 0:2] = Gamma[3, 2:4] = j_bar
            return Gamma, norm(dx_dy)
        elif self.rank == 3:
            i_bar = coords[1, :] - coords[0, :]  # [dx, dy, dz]
            j_bar = np.array([0, 1, 0])  # Assume j-bar points upwards
            k_bar = np.cross(i_bar, j_bar)
            Gamma[np.ix_([0,1,2],[0,1,2])] = gram_schmidt(i_bar, j_bar, k_bar)
            Gamma[np.ix_([3,4,5],[3,4,5])] = Gamma[np.ix_([0, 1, 2], [0, 1, 2])]
            return Gamma, norm(i_bar)
