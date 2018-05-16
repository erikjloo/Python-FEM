# Import Standard Libraries
import re
import scipy as np

#===========================================================================
#   LoadTable
#===========================================================================


class LoadTable(object):
    """ LoadTable """

    # Public:

    #-----------------------------------------------------------------------
    #   constructor
    #-----------------------------------------------------------------------

    def __init__(self):
        self.ndof = 0
        self.loads = np.empty(0)
        self.loads[:] = np.nan

    #-----------------------------------------------------------------------
    #   initialize
    #-----------------------------------------------------------------------

    def initialize(self, props, mesh):
        """ Input:  props = properties """

        self.ndof = mesh.dofCount()
        self.loads = np.empty(self.ndof)
        self.loads[:] = 0

        props = props.getProps("input.loads")
        type = props.get("type")
        path = props.get("file")
        if type == "Input":
            self.readXML(path, mesh)

    #-----------------------------------------------------------------------
    #   readXML
    #-----------------------------------------------------------------------

    def readXML(self, path, mesh):
        """ Input: path = path_to_file """
        with open(path, 'r') as file:

            flag_c = False

            for line in file:
                if line.startswith("<Loads>"):
                    flag_c = True
                elif line.startswith("</Loads>"):
                    flag_c = False

                if flag_c is True and not line.startswith("<Loads>"):
                    dof = re.findall(r"[a-zA-Z]+", line)[0]
                    [node, rval] = re.findall(r"[-+]?\d*\.\d+|\d+", line)
                    # print(dof, "[", node, "] = ", rval)
                    idof = mesh.getDofIndex(int(node), dof)
                    self.addLoad(idof, float(rval))

    #-----------------------------------------------------------------------
    #   addLoad
    #-----------------------------------------------------------------------

    def addLoad(self, idof, rval=0.0):
        self.loads[idof] = rval

    #-----------------------------------------------------------------------
    #   addLoadTable
    #-----------------------------------------------------------------------

    def addLoads(self, idofs, rval=0.0):
        if isinstance(idofs, (list, tuple, range, np.ndarray)):
            for idof in idofs:
                self.addLoad(idof, rval)
        else:
            self.addLoad(idofs, rval)

    def getLoads(self):
        return self.loads

    def dofCount(self):
        return self.ndof
