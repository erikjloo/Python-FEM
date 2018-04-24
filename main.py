# Import Local Libraries
from algebra import MatrixBuilder

from properties import Properties
from factory import ModelFactory


from nonlin import multistep


file = "Examples/square.pro"

props = Properties()
props.parseFile(file)
props.print()

model = ModelFactory(props)
models = model.createModel()



# multistep(model)


# output:
# u, strain, sigma, f_int, K


