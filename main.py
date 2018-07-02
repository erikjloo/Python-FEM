# Import Standard Libraries
import time
import scipy as np

# Import Local Libraries
from properties import Properties
from modules import ChainModule 
from modules import InputModule
from modules import InitModule
from modules import LinSolveModule
from modules import Execute
np.set_printoptions(precision=4)

file = input("Properties file: ")

start = time.time()

# Props & Globdat
conf = Properties()
props = Properties()
globdat = Properties()
props.parseFile(file)

# Chain module
module = ChainModule()

# Create mesh
module.pushBack(InputModule("input"))

# Create model, cons, mbuilder and vectors
module.pushBack(InitModule())

# Linear analysis
module.pushBack(LinSolveModule("linsolve"))

# Execute
Execute(module, conf, props, globdat)

stop = time.time()
print("Elapsed time is ",stop-start)

# globdat.mesh.updateGeometry(disp)
mesh = globdat.get("mesh")
mesh.plotDeformed(globdat.get("solu"), 1)

print(globdat.get("solu"))
