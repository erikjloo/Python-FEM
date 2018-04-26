# Import Standard Libraries
from abc import ABCMeta, abstractclassmethod


class Module(metaclass=ABCMeta):
    """ Abstract Module Class 
    
    Virtual static methods:
        init(name, props, mesh)
    """
    
    @abstractclassmethod
    def init(cls, name, props, mesh):
        # No need to implement
        pass

    @abstractclassmethod
    def run(cls, mesh):
        # No need to implement
        pass

    @abstractclassmethod
    def shutdown(cls, mesh):
        # No need to implement
        pass

class InitModule(Module):

    @classmethod
    def __init__(cls, name, props, mesh):
        cls.__makeModel()
        cls.__makeConstraints()
        cls.__makeVectors()

    @classmethod
    def __makeModel(cls):
        pass

    @classmethod
    def __makeConstraints(cls):
        pass
    
    @classmethod
    def __makeVectors(cls):
        pass


