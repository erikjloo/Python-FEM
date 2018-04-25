# Import Local Libraries

from properties import Properties
from models import ModelFactory
from algebra import MatrixBuilder
from nonlin import multistep


# Initialization
file = "Examples/rve.pro"

props = Properties()
props.parseFile(file)
props.print()

model = ModelFactory(props)
models = model.createModel()



# Solver



# Output



# multistep(model)


# output:
# u, strain, sigma, f_int, K


