# Import Standard Libraries
import scipy as np

# Import Local Libraries
from constraints import Constraints
from properties import Properties
from loadTable import LoadTable
from algebra import MatrixBuilder
# from models import ModelFactory
from mesh import Mesh

class GlobalData(Properties):
    """ Global data 

    Instance Members:
        i = load step number
        ndof = number of degrees of freedom
        mesh = nodeset, elementset and dofspace
        fint = vector of internal forces
        fext = vector of external forces
        disp = solution vector
        load = load table
        cons = constraints
        mbuild = matrix builder
    
    Public Methods:
        makeMesh(props)
        makeModel(props)
        makeLoadTable(props)
        makeConstraints(props)
        makeMatrixBuilder()
        makeVectors()
    """

    def __init__(self, props):
        self.properties = {}
        self.i = 0
        self.ndof = 0
    
