#  Import Standard libraties
import scipy as np

#  Import Local Libraries
from abc import ABCMeta, abstractmethod
from constraints import Constraints
from algebra import MatrixBuilder
from models import ModelFactory
from solvers import Solver
from mesh import Mesh


#===========================================================================
#   Module
#===========================================================================


class Module(metaclass=ABCMeta):
    """ Abstract Module Class 
    
    Pure Virtual Methods:
        init(name, props, globdat)
        run(globdat)
        shutdown(globdat)
    """

    def __init__(self, name=None):
        if name: self.name = name

    @abstractmethod
    def init(self, props, globdat): pass

    @abstractmethod
    def run(self, globdat): pass

    @abstractmethod
    def shutdown(self, globdat): pass


#===========================================================================
#   InputModule
#===========================================================================


class InputModule(Module):
    """ Reads props and initializes the mesh
        
    Public Methods:
        init(name, props, globdat)
        run(globdat)
        shutdown(globdat)
    """

    def init(self, props, globdat):
        props = props.getProps(self.name)
        globdat.makeMesh(props, globdat)

    def run(self, globdat):
        pass

    def shutdown(self, globdat):
        pass

    def __makeMesh(self, props, globdat):
        props = props.getProps("mesh")
        globdat.mesh.initialize(props)

#===========================================================================
#   InitModule
#===========================================================================


class InitModule(Module):
    """ Initializes the model, constraints, matrix builder and vectors
            
    Public Methods:
        init(name, props, globdat)
        run(globdat)
        shutdown(globdat)
    """

    def init(self, props, globdat):
        globdat.makeModel(props, globdat.mesh)
        globdat.makeConstraints(props)
        globdat.makeMatrixBuilder()
        globdat.makeVectors()

    def run(self, globdat):
        pass

    def shutdown(self, globdat):
        pass


#===========================================================================
#   LinSolveModule
#===========================================================================


class LinSolveModule(Module):
    """ Runs a linear analysis
            
    Public Methods:
        init(name, props, globdat)
        run(globdat)
        shutdown(globdat)
    """

    def init(self, props, globdat):
        pass

    def run(self, globdat):

        # Initial stiffness matrix
        hbw = globdat.model.get_Matrix_0(
            globdat.mbuild, globdat.fint, globdat.disp, globdat.mesh)
        print("The half-band-width is", hbw)

        # Solve
        solver = Solver("numpy", globdat.cons)
        K = globdat.mbuild.getDenseMatrix()
        solver.solve(K, globdat.disp, globdat.fext, hbw)
        
    def shutdown(self, globdat):
        pass


#===========================================================================
#   NonlinModule
#===========================================================================
