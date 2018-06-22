# Import Standard Libraries
import scipy as np
import time

# Import Local Libraries
from properties import Properties
from globalData import GlobalData
from modules import ChainModule
from modules import InputModule
from modules import InitModule
from modules import NonlinModule
from modules import SampleModule
from modules import ControlModule
from modules import Execute
np.set_printoptions(precision=4)

start = time.time()

file = input("Properties file: ")

# Props
props = Properties()
conf = Properties()
props.parseFile(file)
globdat = GlobalData(props)

# Chain module
module = ChainModule()
module.pushBack(InputModule("input"))
module.pushBack(InitModule())
module.pushBack(NonlinModule("nonlin"))
module.pushBack(SampleModule("sample"))
module.pushBack(ControlModule("control"))

# Execute
Execute(module, conf, props, globdat)

stop = time.time()
print("Elapsed time is ", stop-start)
globdat.mesh.plotDeformed(globdat.disp,1)
