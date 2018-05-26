# Import Standard Libraries
import scipy as np

# Import Local Libraries
from constraints import Constraints
from loadTable import LoadTable
from algebra import MatrixBuilder
from models import ModelFactory
from mesh import Mesh

class GlobalData(object):
    """ Global data 

    Instance Members:
        ndof = number of degrees of freedom
        mesh = nodeset, elementset and dofspace
        fint = vector of internal forces
        fext = vector of external forces
        disp = vector of displacements (solution)
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
        self.ndof = 0
        self.mesh = Mesh()
        self.fext = np.zeros(0)
        self.fint = np.zeros(0)
        self.disp = np.zeros(0)
        self.load = LoadTable()
        self.cons = Constraints()
        self.mbuild = MatrixBuilder(0)
    
    def makeMesh(self, props):
        self.mesh.initialize(props)

    def makeModel(self, props):
        self.model = ModelFactory("model", props, self.mesh)
        self.ndof = self.mesh.dofCount()

    def makeLoadTable(self, props):
        self.load.initialize(props, self.mesh)

    def makeConstraints(self, props):
        self.cons.initialize(props, self.mesh)

    def makeMatrixBuilder(self):
        self.mbuild.resize(self.ndof)

    def makeVectors(self):
        self.fext.resize(self.ndof)
        self.fint.resize(self.ndof)
        self.disp.resize(self.ndof)
