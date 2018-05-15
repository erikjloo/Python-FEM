#  Import Standard libraties
import scipy as np

#  Import Local Libraries
from abc import ABCMeta, abstractmethod
from constraints import Constraints
from algebra import MatrixBuilder
from models import ModelFactory
from mesh import Mesh


#===========================================================================
#   Module
#===========================================================================


class Module(metaclass=ABCMeta):
    """ Abstract Module Class 
    
    Pure Virtual Methods:
        init(name, props, mesh)
        run(mesh)
        shutdown(mesh)
    """

    @abstractmethod
    def init(self, props, globdat): pass

    @abstractmethod
    def run(self, mesh): pass

    @abstractmethod
    def shutdown(self, mesh): pass


#===========================================================================
#   InputModule
#===========================================================================


class InputModule(Module):
    """ Reads props and initializes the mesh """

    def __init__(self, name):
        self.name = name

    def init(self, props, globdat):
        props = props.getProps(self.name)
        self.__makeMesh(props, globdat)

    def run(self, mesh):
        pass

    def shutdown(self, mesh):
        pass

    def __makeMesh(self, props, globdat):
        props = props.getProps("mesh")
        globdat.mesh = Mesh()
        globdat.mesh.initialize(props)

#===========================================================================
#   InitModule
#===========================================================================


class InitModule(Module):
    """ Initializes the model, constraints, matrix builder and vectors """
    
    def init(self, props, globdat):
        self.__makeModel(props, globdat)
        self.__makeConstraints(props, globdat)
        self.__makeMatrixBuilder(globdat)
        self.__makeVectors(globdat)

    def run(self, mesh):
        pass

    def shutdown(self, mesh):
        pass

    def __makeModel(self, props, globdat):
        globdat.model = ModelFactory("model", props, globdat.mesh)

    def __makeConstraints(self, props, globdat):
        ndof = globdat.mesh.dofCount()
        globdat.cons = Constraints(ndof)
        globdat.cons.initialize(props, globdat.mesh)

    def __makeMatrixBuilder(self, globdat):
        ndof = globdat.mesh.dofCount()
        globdat.mbuild = MatrixBuilder(ndof)

    def __makeVectors(self, globdat):
        ndof = globdat.mesh.dofCount()
        globdat.fint = np.zeros(ndof)
        globdat.fext = np.zeros(ndof)
        globdat.disp = np.zeros(ndof)


#===========================================================================
#   LinSolveModule
#===========================================================================


class LinSolveModule(Module):

    def init(self, props, mesh):
        pass

    def run(self, mesh):
        pass

    def shutdown(self, mesh):
        pass


#===========================================================================
#   NonlinModule
#===========================================================================
