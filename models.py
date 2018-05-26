# Import Standard Libraries
from abc import ABCMeta, abstractmethod


#===========================================================================
#   ModelFactory
#===========================================================================


def ModelFactory(name, props, mesh):
    """ Input:  props = properties
                mesh = mesh
        Output: model """
    props = props.getProps(name)
    type = props.get("type")

    if type == "Matrix":
        # Creates the root
        print("Creating a matrix model named", name)
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
    
    Pure Virtual Methods:
        Model(name, props, mesh)
        get_Matrix_0(mbuild, fint, disp, mesh)
        get_Ext_Vector(fext, mesh)
        get_Int_Vector(fint, disp, mesh)
        get_Constraints(cons, mesh)
        takeAction(action, mesh)
    """

    @abstractmethod
    def __init__(self, name, props, mesh):
        raise NotImplementedError()

    @abstractmethod
    def get_Matrix_0(self, mbuild, fint, disp, mesh):
        raise NotImplementedError()

    @abstractmethod
    def get_Ext_Vector(self, fext, mesh):
        raise NotImplementedError()

    @abstractmethod
    def get_Constraints(self, cons, mesh):
        raise NotImplementedError()

    @abstractmethod
    def takeAction(self, action, mesh):
        raise NotImplementedError()


#===========================================================================
#   MultiModel
#===========================================================================


class MultiModel(Model):
    """ A node in the model tree

    Instance Members:
        name = model name
        models = children

    Public Methods:
        MultiModel(name, props, mesh)
        get_Matrix_0(mbuild, fint, disp, mesh)
        get_Ext_Vector(fext, mesh)
        get_Int_Vector(fint, disp, mesh)
        get_Constraints(cons, mesh)
        takeAction(action, mesh)
    """

    def __init__(self, name, props, mesh):
        """ Creates a node and its children """
        self.name = name

        # Create children
        self.models = []
        for name in props.get("models"):
            model = ModelFactory(name, props, mesh)
            self.models.append(model)

    def get_Matrix_0(self, mbuild, fint, disp, mesh):
        hbw = []
        for model in self.models:
            hbw.append(model.get_Matrix_0(mbuild, fint, disp, mesh))
        return max(hbw)

    def get_Int_Vector(self, fint, disp, mesh):
        for model in self.models:
            model.get_Int_Vector(fint, disp, mesh)

    def get_Ext_Vector(self, fext, mesh):
        for model in self.models:
            model.get_Ext_Vector(fext, mesh)

    def get_Constraints(self, cons, mesh):
        for model in self.models:
            model.get_Constraints(cons, mesh)

    def takeAction(self, action, mesh):
        for model in self.models:
            model.takeAction(action, mesh)

#===========================================================================
#   MatrixModel
#===========================================================================


class MatrixModel(Model):
    """ The root of the model tree 

    Instance Members:
        name = model name
        model = child model

    Public Methods:
        MatrixModel(name, props, mesh)
        get_Matrix_0(mbuild, fint, disp, mesh)
        get_Ext_Vector(fext, mesh)
        get_Int_Vector(fint, disp, mesh)
        get_Constraints(cons, mesh)
        takeAction(action, mesh)
    """

    def __init__(self, name, props, mesh):
        """ Creates a node and its child """
        self.name = name

        # Create child
        self.model = ModelFactory("model", props, mesh)

    def get_Matrix_0(self, mbuild, fint, disp, mesh):
        return self.model.get_Matrix_0(mbuild, fint, disp, mesh)

    def get_Int_Vector(self, fint, disp, mesh):
        self.model.get_Int_Vector(fint, disp, mesh)

    def get_Ext_Vector(self, fext, mesh):
        self.model.get_Ext_Vector(fext, mesh)

    def get_Constraints(self, cons, mesh):
        self.model.get_Constraints(cons, mesh)

    def takeAction(self, action, mesh):
        self.model.takeAction(action, mesh)

# PointLoadModel

# ConstraintsModel

#===========================================================================
#   LoadScaleModel
#===========================================================================


class LoadScaleModel(Model):
    """ The root of the model tree 

    Instance Members:
        name = model name
        model = child model

    Public Methods:
        MatrixModel(name, props, mesh)
        get_Matrix_0(mbuild, fint, disp, mesh)
        get_Ext_Vector(fext, mesh)
        get_Int_Vector(fint, disp, mesh)
        get_Constraints(cons, mesh)
        takeAction(action, mesh)
    """

    def __init__(self, name, props, mesh):
        """ Creates a node and its child """
        self.name = name

        # Create child
        self.model = ModelFactory("model", props, mesh)

    def get_Matrix_0(self, mbuild, fint, disp, mesh):
        return self.model.get_Matrix_0(mbuild, fint, disp, mesh)

    def get_Int_Vector(self, fint, disp, mesh):
        self.model.get_Int_Vector(fint, disp, mesh)

    def get_Ext_Vector(self, fext, mesh):
        self.model.get_Ext_Vector(fext, mesh)

    def get_Constraints(self, cons, mesh):
        self.model.get_Constraints(cons, mesh)

    def takeAction(self, action, mesh):
        self.model.takeAction(action, mesh)


