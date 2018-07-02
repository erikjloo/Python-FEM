# Import Standard Libraries
import time
import scipy as np

# Import Local Libraries
from properties import Properties
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

# Props & Globdat     
conf = Properties()
props = Properties()
globdat = Properties()
props.parseFile(file)

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

mesh = globdat.get("mesh")
mesh.plotDeformed(globdat.get("solu"),1)
