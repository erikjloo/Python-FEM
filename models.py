# Import Standard Libraries
from abc import ABCMeta, abstractmethod


#===========================================================================
#   ModelFactory
#===========================================================================


def ModelFactory(name, props, mesh):

    props = props.getProps(name)
    type = props.get("type")

    if type == "Multi":
        # Creates a node
        print("Creating a multi model named", name)
        return MultiModel(name, props, mesh)

    elif type == "Solid":
        # Creates a leaf
        from solidModel import SolidModel
        print("Creating a solid model named", name)
        return SolidModel(name, props, mesh)

    elif type == "Periodic":
        # Creates a leaf
        from PBCmodel import PBCmodel
        print("Creating a periodic model named", name)
        return PBCmodel(name, props, mesh)


#===========================================================================
#   BaseClass
#===========================================================================


class Model(metaclass=ABCMeta):
    """ Abstract Model Class
    
    """

    @abstractmethod
    def __init__(self, name, props, mesh):
        pass

    @abstractmethod
    def get_Matrix_0(self, mesh, mbuild, f_int):
        pass

    @abstractmethod
    def get_Ext_vector(self, f_ext):
        pass

    @abstractmethod
    def get_Constraints(self, mesh, constraints):
        pass

#===========================================================================
#   MultiModel
#===========================================================================


class MultiModel(object):
    """ A node in the model tree """

    def __init__(self, name, props, mesh):
        """ Creates a node """
        self.name = name

        # Create children
        self.models = []
        for name in props.get("models"):
            model = ModelFactory(name, props, mesh)
            self.models.append(model)

    def initialize(self, props, mesh):
        for model in self.models:
            model.initialize(props, mesh)


#===========================================================================
#   MatrixModel
#===========================================================================


class MatrixModel(object):
    """ The root of the model tree """

    def __init__(self, name, props, mesh):
        """ Creates a node """
        self.name = name

        # Create child model
        name = props.get("model")
        self.model = ModelFactory(name, props, mesh)


#===========================================================================
#   Example
#===========================================================================


if __name__ == "__main__":

    from properties import Properties
    from mesh import Mesh

    file = "Examples/rve.pro"
    props = Properties()
    props.parseFile(file)

    mesh = Mesh()
    mesh.initialize(props, rank=2)

    model = ModelFactory("model", props, mesh)
