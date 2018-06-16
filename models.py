# Import Standard Libraries
from abc import ABCMeta, abstractmethod


#===========================================================================
#   ModelFactory
#===========================================================================

def ModelFactory(name, props, conf, mesh):
    """ Input:  name = model name
                props = properties
                conf = properties
                mesh = mesh
        Output: model """

    type = props.get("{}.type".format(name))

    if type == "Matrix":
        # Creates the root
        print("Creating a matrix model named", name)
        return MatrixModel(name, props, conf, mesh)

    elif type == "Multi":
        # Creates a node
        print("Creating a multi model named", name)
        return MultiModel(name, props, conf, mesh)

    elif type == "Solid":
        # Creates a leaf
        from solidModel import SolidModel
        print("Creating a solid model named", name)
        return SolidModel(name, props, conf, mesh)

    elif type == "Periodic":
        # Creates a leaf
        from PBCmodel import PBCmodel
        print("Creating a periodic model named", name)
        return PBCmodel(name, props, conf, mesh)
    
    else:
        msg = "{} model not yet implemented".format(name)
        raise NotImplementedError(msg)


#===========================================================================
#   Model
#===========================================================================


class Model(metaclass=ABCMeta):
    """ Abstract Model Class
    
    Pure Virtual Methods:
        Model(name, props, conf, mesh)
        takeAction(action, globdat)
    """

    @abstractmethod
    def __init__(self, name, props, conf, mesh):
        raise NotImplementedError()

    @abstractmethod
    def takeAction(self, action, globdat):
        raise NotImplementedError()


#===========================================================================
#   MatrixModel
#===========================================================================


class MatrixModel(Model):
    """ The root of the model tree 

    Instance Members:
        name = model name
        model = child model

    Public Methods:
        MatrixModel(name, props, conf, mesh)
        takeAction(action, globdat)
    """

    def __init__(self, name, props, conf, mesh):
        """ Creates a node and its child """
        self.name = name
        myProps = props.getProps(name)
        myConf = conf.makeProps(name)

        self.type = myProps.get("type")
        myConf.set("type", self.type)

        # Create child
        self.model = ModelFactory("model", myProps, myConf, mesh)

    def takeAction(self, action, globdat):
        self.model.takeAction(action, globdat)


#===========================================================================
#   MultiModel
#===========================================================================


class MultiModel(Model):
    """ A node in the model tree

    Instance Members:
        name = model name
        models = children

    Public Methods:
        MultiModel(name, props, conf, mesh)
        takeAction(action, globdat)
    """

    def __init__(self, name, props, conf, mesh):
        """ Creates a node and its children """
        self.name = name
        myProps = props.getProps(name)
        myConf = conf.makeProps(name)

        self.type = myProps.get("type")
        sub_models = myProps.get("models")

        myConf.set("type",self.type)
        myConf.set("models",sub_models)

        # Create children
        self.models = []
        for name in sub_models:
            model = ModelFactory(name, myProps, myConf, mesh)
            self.models.append(model)

    def takeAction(self, action, globdat):
        for model in self.models:
            model.takeAction(action, globdat)
            
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
        MatrixModel(name, props, conf, mesh)
        takeAction(action, globdat)
    """

    def __init__(self, name, props, conf, mesh):
        """ Creates a node and its child """
        self.name = name
        myProps = props.getProps(name)
        myConf = conf.makeProps(name)

        # Create child
        self.model = ModelFactory("model", myProps, myConf, mesh)

    def takeAction(self, action, globdat):
        self.model.takeAction(action, globdat)
