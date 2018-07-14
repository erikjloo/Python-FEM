
import os
import logging
from properties import Properties
from modules import ChainModule, InitModule, InputModule, Execute
directory = "_mshes"

# Props & Globdat
props = Properties()
conf = Properties()
globdat = Properties()
props.parseFile("Examples/rve.pro")
coarsen_factors = [0.01592, 0.03, 0.04, 0.05, 0.07, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0]

for file in os.listdir(directory):
    for cf in coarsen_factors:
        filename = os.fsdecode(file)
        pr = props.getProps("input.mesh")
        filename = pr.set("file", "_mshes/"+filename)
        pr = props.getProps("model.pbc")
        pr.set("coarsenFactor",cf)
        print(filename)

        # Chain module
        module = ChainModule()
        module.pushBack(InputModule("input"))
        module.pushBack(InitModule())
        module.init(conf, props, globdat)
