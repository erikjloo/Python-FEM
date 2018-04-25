# Import Standard Libraries
import scipy as np
from abc import ABC, abstractmethod

class Module(ABC):
    """ Abstract Module Class """
    
    @abstractmethod
    def __init__(self, name, props, mesh):
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

class InitModule(Module):

    def __init__(self, name, props, mesh):
        self.__makeModel()
        self.__makeConstraints()
        self.__makeVectors()

    def __makeModel(self):
        pass

    def __makeConstraints(self):
        pass
    
    def __makeVectors(self):
        pass

class InputModule(Module):

    def __init__(self, name, props, globdat):
        myProps = props.findProps(name)
        pass
    
