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
conf = Properties()
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
module.init(props, conf, globdat)

conf.print()

globdat.model.takeAction("PLOT_BOUNDARY",globdat)

module.run(globdat)
module.run(globdat)

module.shutdown(globdat)

stop = time.time()
print("Elapsed time is ", stop-start)
globdat.mesh.plotDeformed(globdat.disp,1)
# globdat.mesh.plotMesh()
