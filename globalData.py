# Import Standard Libraries
import scipy as np

# Import Local Libraries
from constraints import Constraints
from properties import Properties
from loadTable import LoadTable
from algebra import MatrixBuilder
from models import ModelFactory
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
        self.i = 0
        self.ndof = 0
        self.mesh = Mesh()
        self.fext = np.zeros(0)
        self.fint = np.zeros(0)
        self.disp = np.zeros(0)
        self.load = LoadTable()
        self.cons = Constraints()
        self.mbuild = MatrixBuilder(0)
    
    def makeMesh(self, props, conf):
        self.mesh.initialize(props, conf)

    def makeModel(self, props, conf):
        self.model = ModelFactory("model", props, conf, self.mesh)
        self.ndof = self.mesh.dofCount()

    def makeLoadTable(self, props):
        # self.set("load",self.load)
        self.load.initialize(props, self.mesh)

    def makeConstraints(self, props):
        # self.set("cons",self.cons)
        self.cons.initialize(props, self.mesh)

    def makeMatrixBuilder(self):
        self.mbuild.resize(self.ndof)

    def makeVectors(self):
        self.fext.resize(self.ndof)
        self.fint.resize(self.ndof)
        self.disp.resize(self.ndof)
