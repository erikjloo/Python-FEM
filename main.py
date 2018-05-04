# Import Standard Libraries
import scipy as np

# Import Local Libraries
from properties import Properties
from algebra import MatrixBuilder
from modules import InputModule
from models import ModelFactory
from nonlin import multistep
from mesh import Mesh

# Initialization
file = "Examples/rve.pro"
props = Properties()
props.parseFile(file)
mesh = Mesh()

module = InputModule("input")
module.init(props, mesh)

model = ModelFactory("model", props, mesh)
model.takeAction("plot_boundary", mesh)

ndof = mesh.dofCount()
mbuild = MatrixBuilder(ndof)
f_int = np.zeros(ndof)

model.get_Matrix_0(mbuild, f_int, mesh)

# Modules works as a series of commands

def execute(module):
    module.init()
    module.run()
    module.shutdown()

# if __name__ == "__main__":
#     execute(module)
