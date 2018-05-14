#  Import Standard libraties
import scipy as np

#  Import Local Libraries
from abc import ABCMeta, abstractmethod
from constraints import Constraints
from models import ModelFactory
from mesh import Mesh


#===========================================================================
#   Module
#===========================================================================


class Module(metaclass=ABCMeta):
    """ Abstract Module Class 
    
    Instance members:
        name = string for finding runtime parameters

    Public Methods:
        __init__(name)

    Pure Virtual Methods:
        init(name, props, mesh)
        run(mesh)
        shutdown(mesh)
    """
    
    @abstractmethod
    def __init__(self, name):
        pass

    @abstractmethod
    def init(self, props, mesh):
        # No need to implement
        pass

    @abstractmethod
    def run(self, mesh):
        # No need to implement
        pass

    @abstractmethod
    def shutdown(self, mesh):
        # No need to implement
        pass


#===========================================================================
#   InputModule
#===========================================================================


class InputModule(Module):

    def __init__(self, name):
        self.name = name

    def init(self, props):
        props = props.getProps(self.name)
        mesh = self.__makeMesh(props)
        # model = self.__makeModel(props, mesh)
        # cons = self.__makeConstraints(props, mesh)
        # fint, fext, disp = self.__makeVectors(mesh)
        # return mesh, model, cons, fint, fext, disp
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

    # def __makeModel(self, props, mesh):
    #     model = ModelFactory("model", props, mesh)
    #     return model

    # def __makeConstraints(self, props, mesh):
    #     props = props.getProps("input.constraints")
    #     ndof = mesh.dofCount()
    #     cons = Constraints(ndof)
    #     cons.initialize(props, mesh)
    #     return cons

    # def __makeVectors(self, mesh):
    #     ndof = mesh.dofCount()
    #     fint = np.zeros(ndof)
    #     fext = np.zeros(ndof)
    #     disp = np.zeros(ndof)
    #     return fint, fext, disp


#===========================================================================
#   LinSolveModule
#===========================================================================
