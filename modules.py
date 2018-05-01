# Import Standard Libraries
from abc import ABCMeta, abstractmethod


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

    def init(self, props, mesh):
        props = props.getProps(self.name)
        self.__makeMesh(props, mesh)

    def run(self, mesh):
        pass

    def shutdown(self, mesh):
        pass

    def __makeMesh(self, props, mesh):
        props = props.getProps("mesh")
        mesh.initialize(props)


#===========================================================================
#   InitModule
#===========================================================================


class InitModule(Module):

    def __init__(self, name):
        self.name = name

    def init(self, props, mesh):
        props = props.getProps(self.name)
        self.__makeModel(props, mesh)
        self.__makeConstraints(props, mesh)
        self.__makeVectors(props, mesh)

    def run(self, mesh):
        pass

    def shutdown(self, mesh):
        pass

    def __makeModel(self, props, mesh):
        pass

    def __makeConstraints(self, props, mesh):
        pass
    
    def __makeVectors(self, props, mesh):
        pass


#===========================================================================
#   LinSolveModule
#===========================================================================
