# Import Standard Libraries
import time
import logging
import argparse
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

def main():
    logging.basicConfig(level=logging.DEBUG, format='%(message)s')

    parser = argparse.ArgumentParser()
    parser.add_argument("file_path", help="Path to properties file", type = str)
    file = parser.parse_args().file_path
    start = time.time()

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
    logging.info("Elapsed time is {}".format(stop-start))

    mesh = globdat.get("mesh")
    mesh.plotDeformed(globdat.get("solu"),1)

if __name__ == "__main__":
    main()
