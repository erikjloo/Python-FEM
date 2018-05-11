# Import Standard Libraries
import scipy as np

# Import Local Libraries
from properties import Properties
from constraints import Constraints
from algebra import MatrixBuilder
from modules import InputModule
from models import ModelFactory
from nonlin import multistep
from mesh import Mesh

# Initialization
file = "Examples/semicircle.pro"
props = Properties()
props.parseFile(file)
mesh = Mesh()
module = InputModule("input")
module.init(props, mesh)

model = ModelFactory("model", props, mesh)

ndof = mesh.dofCount()
cons = Constraints(ndof)
mbuild = MatrixBuilder(ndof)
fint = np.zeros(ndof)
fext = np.zeros(ndof)
disp = np.zeros(ndof)

idofs = mesh.getDofIndices([1, 4],["u","v","w"])
cons.addConstraints(idofs)
idofs = mesh.getDofIndices([2, 5], ["v", "w"])
cons.addConstraints(idofs)
idofs = mesh.getDofIndices([3, 6], ["u", "w"])

model.get_Matrix_0(mbuild, fint, disp, mesh)

# Modules works as a series of commands

def execute(module):
    module.init()
    module.run()
    module.shutdown()

# if __name__ == "__main__":
#     execute(module)
