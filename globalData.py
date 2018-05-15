# Import Standard Libraries
import scipy as np

# Import Local Libraries
from properties import  Properties
# from constraints import Constraints
# from algebra import MatrixBuilder
# from models import Model
# from mesh import Mesh

class GlobalData(Properties):
    """ Global data """

    def __init__(self):
        self.globdat = {}
        # self.mesh = Mesh()
        # self.fint = np.zeros(0)
        # self.fext = np.zeros(0)
        # self.disp = np.zeros(0)
        # self.cons = Constraints(0)
        # self.mbuild = MatrixBuilder(0)
    
