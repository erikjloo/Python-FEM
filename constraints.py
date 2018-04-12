# Import Standard Libraries
import scipy as np
import warnings
from scipy import ix_

# Import Local Libraries
from dofspace import DofSpace

#===========================================================================
#   Constraints
#===========================================================================


class Constraints(object):


    def __init__(self, dofs):
        if isinstance(dofs, DofSpace):
            self.dofs = dofs