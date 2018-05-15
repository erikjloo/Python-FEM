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

dim = 2
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
if dim == 2:
    idofs = globdat.mesh.getDofIndices(0, ["u", "v"])
    globdat.cons.addConstraints(idofs)
    idof = globdat.mesh.getDofIndex(2, "v")
    globdat.cons.addConstraint(idof)
    idof = globdat.mesh.getDofIndex(1, "u")
    globdat.cons.addConstraint(idof)

    # Loads
    idof = globdat.mesh.getDofIndex(1, "v")
    globdat.fext[idof] = 10

elif dim == 3:
    idofs = globdat.mesh.getDofIndices([0, 3], ["u", "v", "w"])
    globdat.cons.addConstraints(idofs)
    idofs = globdat.mesh.getDofIndices([1, 4], ["v", "w"])
    globdat.cons.addConstraints(idofs)
    idofs = globdat.mesh.getDofIndices([2, 5], ["u", "w"])
    globdat.cons.addConstraints(idofs)

    # Loads
    idofs = globdat.mesh.getDofIndices([2, 5], "v")
    globdat.fext[idofs] = 5

module = LinSolveModule()
module.run(globdat)

stop = time.time()
print("Elapsed time is ",stop-start)

# globdat.mesh.updateGeometry(disp)
globdat.mesh.plotDeformed(globdat.disp, 21.5)