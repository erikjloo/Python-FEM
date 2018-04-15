# Import Standard Libraries
import os
import scipy as np
from configparser import ConfigParser

# Import Local Libraries
from dofspace import DofSpace
from solidModel import SolidModel
from constraints import Constraints
from properties import Properties, ElementType
from algebra import MatrixBuilder


path = '/home/erik/Documents/Python/FEM/rve.msh'


print(("Mesh read with {} nodes and {} elements").format(nnod,nele))

config = ConfigParser()
config.read('square.pro')
ConfigParser
u = np.empty(ndof)

shape = Tri3()


mbuild = MatrixBuilder(ndof)

solid = SolidModel()
solid.readMesh(path)
solid.initialize(2)
rank = 2

K,f = solid.assemble(mbuild)

# output:
# u, strain, sigma, f_int, K


