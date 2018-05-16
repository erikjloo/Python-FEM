# Import Standard Libraries
import matplotlib.pyplot as plt
import scipy as np
import time

# Import Local Libraries
from properties import Properties
from globalData import GlobalData
from modules import InputModule, InitModule, LinSolveModule
from solvers import Solver

np.set_printoptions(precision=4)

start = time.time()

dim = 3
# Initialization
if dim == 2:
    file = "Examples/2D_semicircle.pro"
elif dim == 3:
    file = "Examples/3D_semicircle.pro"

# Props
props = Properties()
props.parseFile(file)
props.print()

# Global Data
globdat = GlobalData(props)

# Create mesh
module = InputModule("input")
module.init(props, globdat)

# Create model, cons, mbuilder and vectors
module = InitModule()
module.init(props, globdat)

# Constraints


module = LinSolveModule()
module.init(props, globdat)
module.run(globdat)

stop = time.time()
print("Elapsed time is ", stop-start)

# globdat.mesh.updateGeometry(disp)
globdat.mesh.plotDeformed(globdat.disp, 21.5)
