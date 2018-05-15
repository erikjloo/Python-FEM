# Import Standard Libraries
import scipy as np

# Import Local Libraries
from properties import  Properties
from constraints import Constraints
from algebra import MatrixBuilder
from models import ModelFactory
from mesh import Mesh

class GlobalData(Properties):
    """ Global data """

    def __init__(self, props):
        self.globdat = props
        self.mesh = Mesh()
        self.fint = np.zeros(0)
        self.fext = np.zeros(0)
        self.disp = np.zeros(0)
        self.cons = Constraints(0)
        self.mbuild = MatrixBuilder(0)
    
    def makeMesh(self, props):
        props = props.getProps("mesh")
        self.mesh.initialize(props)
        
    def makeModel(self, props, mesh):
        self.model = ModelFactory("model", props, mesh)

    def makeConstraints(self, props):
        ndof = self.mesh.dofCount()
        self.cons = Constraints(ndof)
        self.cons.initialize(props, self.mesh)

    def makeMatrixBuilder(self):
        ndof = self.mesh.dofCount()
        self.mbuild = MatrixBuilder(ndof)

    def makeVectors(self):
        ndof = self.mesh.dofCount()
        self.fint = np.zeros(ndof)
        self.fext = np.zeros(ndof)
        self.disp = np.zeros(ndof)
