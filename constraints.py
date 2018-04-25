# Import Standard Libraries
import scipy as np
import warnings
from scipy import ix_


#===========================================================================
#   Constraints
#===========================================================================


class Constraints(object):
    """ Constraints
    
    Static Members:
        __len_error__ = "jdofs and coeff are not the same length!"
    
    """
    __len_error__ = "jdofs and coeff are not the same length!"

    def __init__(self, ndof):
        self.ndof = ndof
        self.conspace = np.nan((ndof,ndof))

    def addConstraint(self, idof, rval=0.0, jdofs=None, coeff=None):
        if jdofs is None:
            self.conspace[idof,0] = rval
        elif len(jdofs) == len(coeff):
            self.conspace[idof,ix_(jdofs)] = coeff
        else:
            warnings.warn(self.__len_error__)
    
    def get_Dirichlet(self):
        return self.conspace

    def dofCount(self):
        return self.ndof

    def get_fdof(self):
        return np.argwhere(np.isnan(self.conspace))
    
    def get_sdof(self):
        return np.argwhere(np.isfinite(self.conspace))
