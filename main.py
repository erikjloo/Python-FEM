# Import Standard Libraries
import scipy as np
import time

# Import Local Libraries
from properties import Properties
from globalData import GlobalData
from modules import ChainModule, InputModule, InitModule, LinSolveModule
# np.set_printoptions(precision=4)

file = input("Properties file: ")

start = time.time()

# Properties
props = Properties()
props.parseFile(file)

# Global Data
globdat = GlobalData(props)

# Chain module
module = ChainModule()

# Create mesh
module.pushBack(InputModule("input"))

# Create model, cons, mbuilder and vectors
module.pushBack(InitModule())

# Linear analysis
module.pushBack(LinSolveModule("linsolve"))

# Execute
module.init(props, globdat)
module.run(globdat)
module.shutdown(globdat)

stop = time.time()
print("Elapsed time is ",stop-start)

# globdat.mesh.updateGeometry(disp)
globdat.mesh.plotDeformed(globdat.disp, 1)

print(globdat.disp)
