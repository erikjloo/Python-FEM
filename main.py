# Import Local Libraries
from algebra import MatrixBuilder
from dofspace import DofSpace
from properties import Properties, ElementType
from solidModel import SolidModel
from constraints import Constraints
from nonlin import multistep


file = "Examples/square.pro"

props = Properties()
props.parseFile(file)
props.print()

userInput = props.getProps("userInput")
models = props.getProps("models")

model = modelFactory(props)

multistep(model)


# output:
# u, strain, sigma, f_int, K


