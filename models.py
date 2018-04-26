# Import Standard Libraries
from abc import ABCMeta, abstractmethod


#===========================================================================
#   ModelFactory
#===========================================================================


def ModelFactory(name, props, mesh):

    props = props.getProps(name)
    type = props.get("type")

    if type == "Matrix":
        # Creates the root
        print("Creating a matrix model names", name)
        return MatrixModel(name, props, mesh)

    elif type == "Multi":
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
#   Model
#===========================================================================


class Model(metaclass=ABCMeta):
    """ Abstract Model Class
    
    Virtual Methods:
        Model(name, props, mesh)
        get_Matrix_0(mbuild, f_int, mesh)
        get_Ext_Vector(f_ext)
        get_Constraints(mesh, constraints)
    """

    @abstractmethod
    def __init__(self, name, props, mesh):
        pass

    @abstractmethod
    def get_Matrix_0(self, mbuild, f_int, mesh):
        pass

    @abstractmethod
    def get_Ext_Vector(self, f_ext):
        pass

    @abstractmethod
    def get_Constraints(self, mesh, constraints):
        pass


#===========================================================================
#   MultiModel
#===========================================================================


class MultiModel(Model):
    """ A node in the model tree

    Instance Members:
        name = model name
        models = children

    Public Methods:
        Model(name, props, mesh)
        get_Matrix_0(mbuild, f_int, mesh)
        get_Ext_Vector(f_ext)
        get_Constraints(mesh, constraints)
    """
    def __init__(self, name, props, mesh):
        """ Creates a node and its children"""
        self.name = name

        # Create children
        self.models = []
        for name in props.get("models"):
            model = ModelFactory(name, props, mesh)
            self.models.append(model)

    def get_Matrix_0(self, mbuild, f_int, mesh):
        for model in self.models:
            model.get_Matrix_0(mbuild, f_int, mesh)

    def get_Ext_Vector(self, f_ext):
        for model in self.models:
            model.get_Ext_Vector(f_ext)

    def get_Constraints(self, mesh, constraints):
        for model in self.models:
            model.get_Constraints(mesh, constraints)


#===========================================================================
#   MatrixModel
#===========================================================================


class MatrixModel(Model):
    """ The root of the model tree 

    Instance Members:
        name = model name
        model = child model

    Public Methods:
        Model(name, props, mesh)
        get_Matrix_0(mbuild, f_int, mesh)
        get_Ext_Vector(f_ext)
        get_Constraints(mesh, constraints)
    """

    def __init__(self, name, props, mesh):
        """ Creates a node and its child """
        self.name = name

        # Create child model
        name = props.get("model")
        self.model = ModelFactory(name, props, mesh)

    def get_Matrix_0(self, mbuild, f_int, mesh):
        self.model.get_Matrix_0(mbuild, f_int, mesh)

    def get_Ext_Vector(self, f_ext):
        self.model.get_Ext_Vector(f_ext)

    def get_Constraints(self, mesh, constraints):
        self.model.get_Constraints(mesh, constraints)

# LoadScaleModel

# PointLoadModel

# ConstraintsModel

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
    mesh.initialize(props.getProps("input.mesh"))

    model = ModelFactory("model", props, mesh)
