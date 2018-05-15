# Import Standard Libraries
import scipy as np
import warnings

#===========================================================================
#   Constraints
#===========================================================================


class Constraints(object):
    """ Constraints """

    def __init__(self, ndof):
        self.ndof = ndof
        self.conspace = np.empty(ndof)
        self.conspace[:] = np.nan

    def initialize(self, props, mesh):
        pass

    def addConstraint(self, idof, rval=0.0):
        self.conspace[idof] = rval
    
    def addConstraints(self, idofs, rval=0.0):
        if isinstance(idofs,(list,tuple,range,np.ndarray)):
            for idof in idofs:
                self.addConstraint(idof, rval)
        else:
            self.addConstraint(idofs, rval)

    def getConspace(self):
        return self.conspace

    def dofCount(self):
        return self.ndof

    def get_fdof(self):
        return np.argwhere(np.isnan(self.conspace)).transpose()[0]
    
    def get_sdof(self):
        return np.argwhere(np.isfinite(self.conspace)).transpose()[0]
