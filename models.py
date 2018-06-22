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

    def __init__(self, name, props, conf, mesh):
        self.name = name

    def __del__(self):
        print("Cleaning {} model".format(self.name))

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

        myConf.set("type", self.type)
        myConf.set("models", sub_models)

        # Create children
        self.models = []
        for name in sub_models:
            model = ModelFactory(name, myProps, myConf, mesh)
            self.models.append(model)

    def takeAction(self, action, globdat):
        for model in self.models:
            model.takeAction(action, globdat)


#===========================================================================
#   PointLoadModel
#===========================================================================


class PointLoadModel(Model):
    """ Assigns forces to the external force vector

    Instance Members:
        name = model name

    Public Methods:
        PointLoadModel(name, props, conf, mesh)
        takeAction(action, globdat)
    """

    __key__ = "PointLoadModel: loadTable not specified!"

    def __init__(self, name, props, conf, mesh):
        self.name = name
        myProps = props.getProps(name)
        myConf = conf.makeProps(name)

        self.loadTable = myProps.get("loadTable")
        myConf.set("loadTable", self.loadTable)

        if self.loadTable is None:
            raise KeyError(self.__key__)
            
    def takeAction(self, action, globdat):
        if action == "GET_EXT_VECTOR":
            loadTable = globdat.get(self.loadTable)
            loadScale = globdat.get("loadScale")
            globdat.fext += loadScale*loadTable.getLoads()
            return True
        else:
            return False


#===========================================================================
#   ConstraintsModel
#===========================================================================


class ConstraintsModel(Model):
    """ Assigns constraints

    Instance Members:
        name = model name

    Public Methods:
        ConstraintsModel(name, props, conf, mesh)
        takeAction(action, globdat)
    """

    __key__ = "ConstraintsModel: constraints not specified!"

    def __init__(self, name, props, conf, mesh):
        self.name = name
        myProps = props.getProps(name)
        myConf = conf.makeProps(name)

        self.constraints = myProps.get("constraints")
        myConf.set("constraints", self.constraints)

        if self.constraints is None:
            raise KeyError(self.__key__)

    def takeAction(self, action, globdat):
        if action == "GET_CONSTRAINTS":
            constraints = globdat.get(self.constraints)
            loadScale = globdat.get("loadScale")

            return True
        else:
            return False

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
