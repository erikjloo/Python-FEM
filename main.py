# Import Local Libraries

from properties import Properties
from algebra import MatrixBuilder
from models import ModelFactory
from nonlin import multistep
from mesh import Mesh

# Initialization
file = "Examples/rve.pro"
props = Properties()
props.parseFile(file)

mesh = Mesh()
mesh.initialize(props, rank=2)

model = ModelFactory("model", props, mesh)

# Solver



# Output



# multistep(model)


# output:
# u, strain, sigma, f_int, K


