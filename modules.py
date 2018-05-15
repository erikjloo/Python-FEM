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
    def init(self, props, mesh): pass

    @abstractmethod
    def run(self, mesh): pass

    @abstractmethod
    def shutdown(self, mesh): pass


#===========================================================================
#   InputModule
#===========================================================================


class InputModule(Module):

    def __init__(self, name):
        self.name = name

    def init(self, props):
        props = props.getProps(self.name)
        mesh = self.__makeMesh(props)
        return mesh

    def run(self, mesh):
        pass

    def shutdown(self, mesh):
        pass

    def __makeMesh(self, props):
        mesh = Mesh()
        props = props.getProps("mesh")
        mesh.initialize(props)
        return mesh

#===========================================================================
#   InitModule
#===========================================================================


class InitModule(Module):

    def init(self, props, mesh):
        model = self.__makeModel(props, mesh)
        cons = self.__makeConstraints(props, mesh)
        mbuild = self.__makeMatrixBuilder(mesh)
        [fint, fext, disp] = self.__makeVectors(mesh)
        return [model, cons, mbuild, fint, fext, disp]

    def run(self, mesh):
        pass

    def shutdown(self, mesh):
        pass

    def __makeModel(self, props, mesh):
        model = ModelFactory("model", props, mesh)
        return model

    def __makeConstraints(self, props, mesh):
        ndof = mesh.dofCount()
        cons = Constraints(ndof)
        cons.initialize(props, mesh)
        return cons

    def __makeMatrixBuilder(self, mesh):
        ndof = mesh.dofCount()
        mbuild = MatrixBuilder(ndof)
        return mbuild

    def __makeVectors(self, mesh):
        ndof = mesh.dofCount()
        fint = np.zeros(ndof)
        fext = np.zeros(ndof)
        disp = np.zeros(ndof)
        return [fint, fext, disp]


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
