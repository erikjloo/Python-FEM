# Import Standard Libraries
import re
import scipy as np

#===========================================================================
#   Constraints
#===========================================================================


class Constraints(object):
    """ Constraints """

    # Public:

    #-----------------------------------------------------------------------
    #   constructor
    #-----------------------------------------------------------------------

    def __init__(self):
        self.ndof = 0
        self.conspace = np.empty(0)
        self.conspace[:] = np.nan

    #-----------------------------------------------------------------------
    #   initialize
    #-----------------------------------------------------------------------

    def initialize(self, props, mesh):
        """ Input:  props = properties """

        self.ndof = mesh.dofCount()
        self.conspace = np.empty(self.ndof)
        self.conspace[:] = np.nan
        
        try:
            props = props.getProps("input.constraints")
            if props.get("type") == "Input":
                path = props.get("file")
                self.readXML(path, mesh)
                print(path,"file read")
        except:
            KeyError(" No constraints provided ")

    #-----------------------------------------------------------------------
    #   readXML
    #-----------------------------------------------------------------------

    def readXML(self, path, mesh):
        """ Input: path = path_to_file """
        with open(path, 'r') as file:
            
            flag_c = False

            for line in file:
                if line.startswith("<Constraints>"):
                    flag_c = True
                elif line.startswith("</Constraints>"):
                    flag_c = False
                
                if flag_c is True and not line.startswith("<Constraints>"):
                    dof = re.findall(r"[a-zA-Z]+", line)[0]
                    [node, rval] = re.findall(r"[-+]?\d*\.\d+|\d+", line)
                    # print(dof, "[", node, "] = ", rval)
                    idof = mesh.getDofIndex(int(node), dof)
                    self.addConstraint(idof, float(rval))

    #-----------------------------------------------------------------------
    #   addConstraint
    #-----------------------------------------------------------------------

    def addConstraint(self, idof, rval=0.0):
        self.conspace[idof] = rval

    #-----------------------------------------------------------------------
    #   addConstraints
    #-----------------------------------------------------------------------

    def addConstraints(self, idofs, rval=0.0):
        if isinstance(idofs,(list,tuple,range,np.ndarray)):
            for idof in idofs:
                self.addConstraint(idof, rval)
        else:
            self.addConstraint(idofs, rval)

    #-----------------------------------------------------------------------
    #   getConspace
    #-----------------------------------------------------------------------

    def getConspace(self):
        return self.conspace

    #-----------------------------------------------------------------------
    #   dofCount
    #-----------------------------------------------------------------------

    def dofCount(self):
        return self.ndof

    #-----------------------------------------------------------------------
    #   get_fdof
    #-----------------------------------------------------------------------

    def get_fdof(self):
        return np.argwhere(np.isnan(self.conspace)).transpose()[0]

    #-----------------------------------------------------------------------
    #   get_sdof
    #-----------------------------------------------------------------------

    def get_sdof(self):
        return np.argwhere(np.isfinite(self.conspace)).transpose()[0]
