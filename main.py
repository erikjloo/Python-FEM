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
mesh.initialize(props)

model = ModelFactory("model", props, mesh)

# Modules works as a series of commands


def execute(module):
    module.init()
    module.run()
    module.shutdown()

if __name__ == "__main__":
    execute(module)
# Solver



# Output



# multistep(model)


# output:
# u, strain, sigma, f_int, K


