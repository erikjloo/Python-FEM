
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
coarsen_factors = [0.01592, 0.025, 0.04, 0.063, 0.1, 0.16, 0.25, 0.4, 0.63, 1.0]

with open('Examples/trCount.dat', 'w') as f:
    f.write(" {} \n".format(coarsen_factors))

for file in os.listdir(directory):
    filename = os.fsdecode(file)
    pr = props.getProps("input.mesh")
    filename = pr.set("file", "_mshes/"+filename)
    print(filename)
    with open('Examples/trCount.dat', 'a') as f:
        f.write(filename)

    for cf in coarsen_factors:
        pr = props.getProps("model.pbc")
        pr.set("coarsenFactor",cf)
        module = ChainModule()
        module.pushBack(InputModule("input"))
        module.pushBack(InitModule())
        module.init(conf, props, globdat)
        # Modify PBCmodel so it outputs len(trNodes) to file

    with open('Examples/trCount.dat', 'a') as f:
        f.write(" \n ")
