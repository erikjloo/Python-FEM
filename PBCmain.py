# Import Standard Libraries
import scipy as np
import time

# Import Local Libraries
from properties import Properties
from globalData import GlobalData
from modules import ChainModule, InputModule, InitModule, NonlinModule

np.set_printoptions(precision=4)

start = time.time()

file = "Examples/square.pro"

# Props
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
module.pushBack(NonlinModule("nonlin"))

# Execute
module.init(props, globdat)
globdat.model.takeAction("plot_boundary",globdat.mesh)

for _ in range(1):
    module.run(globdat)

module.shutdown(globdat)

stop = time.time()
print("Elapsed time is ", stop-start)
globdat.mesh.plotDeformed(globdat.disp,5)
# globdat.mesh.plotMesh()
