# Import Standard Libraries
from abc import ABCMeta, abstractmethod


#===========================================================================
#   ModelFactory
#===========================================================================


def ModelFactory(name, conf, props, mesh):
    """ Input:  name = model name
                conf = output properties
                props = input properties
                mesh = mesh
        Output: model """

    message = "Creating a {} model named {}"
    type = props.get("{}.type".format(name))

    if type == "Matrix":
        # Creates the root
        print(message.format(type, name))
        return MatrixModel(name, conf, props, mesh)

    elif type == "Multi":
        # Creates a node
        print(message.format(type, name))
        return MultiModel(name, conf, props, mesh)

    elif type == "Solid":
        # Creates a leaf
        from solidModel import SolidModel
        print(message.format(type, name))
        return SolidModel(name, conf, props, mesh)

    elif type == "Periodic":
        # Creates a leaf
        from PBCmodel import PBCmodel
        print(message.format(type, name))
        return PBCmodel(name, conf, props, mesh)

    elif type == "PointLoad":
        # Creates a leaf
        print(message.format(type, name))
        return PointLoadModel(name, conf, props, mesh)

    elif type == "Constraints":
        # Creates a leaf
        print(message.format(type, name))
        return ConstraintsModel(name, conf, props, mesh)

    elif type == "LoadScale":
        # Creates a node
        print(message.format(type, name))
        return LoadScaleModel(name, conf, props, mesh)

    else:
        msg = "{} model not implemented".format(name)
        raise KeyError(msg)


#===========================================================================
#   Model
#===========================================================================


class Model(metaclass=ABCMeta):
    """ Abstract Model Class
    
    Pure Virtual Methods:
        Model(name, conf, props, mesh)
        takeAction(action, globdat)
    """

    __key__ = "{}: {} not specified!"

    def __init__(self, name, conf, props, mesh):
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
        type = model type ("Matrix")
        model = child model

    Public Methods:
        MatrixModel(name, conf, props, mesh)
        takeAction(action, globdat)
    """

    def __init__(self, name, conf, props, mesh):
        """ Creates a node and its child """
        self.name = name
        myProps = props.getProps(name)
        myConf = conf.makeProps(name)

        self.type = myProps.get("type", "Matrix")
        myConf.set("type", self.type)

        # Create child
        self.model = ModelFactory("model", myConf, myProps, mesh)

    def takeAction(self, action, globdat):
        self.model.takeAction(action, globdat)


#===========================================================================
#   MultiModel
#===========================================================================


class MultiModel(Model):
    """ A node in the model tree

    Instance Members:
        name = model name
        type = model type ("Multi")
        models = children

    Public Methods:
        MultiModel(name, conf, props, mesh)
        takeAction(action, globdat)
    """

    def __init__(self, name, conf, props, mesh):
        """ Creates a node and its children """
        self.name = name
        myProps = props.getProps(name)
        myConf = conf.makeProps(name)

        self.type = myProps.get("type", "Multi")
        sub_models = myProps.get("models")

        myConf.set("type", self.type)
        myConf.set("models", sub_models)

        if sub_models is None:
            raise KeyError(self.__key__.format(self.type, "models"))

        # Create children
        self.models = []
        for name in sub_models:
            model = ModelFactory(name, myConf, myProps, mesh)
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
        type = model type("PointLoad")
        loadTable = name of LoadTable object

    Public Methods:
        PointLoadModel(name, conf, props, mesh)
        takeAction(action, globdat)
    """

    def __init__(self, name, conf, props, mesh):
        self.name = name
        myProps = props.getProps(name)
        myConf = conf.makeProps(name)

        self.type = myProps.get("type", "PointLoad")
        self.loadTable = myProps.get("loadTable")

        myConf.set("type", self.type)
        myConf.set("loadTable", self.loadTable)

        if self.loadTable is None:
            raise KeyError(self.__key__.format(self.type, "loadTable"))

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
        type = model type ("Constraints")
        conTable = name of Constraints object

    Public Methods:
        ConstraintsModel(name, conf, props, mesh)
        takeAction(action, globdat)
    """

    __key__ = "ConstraintsModel: constraints not specified!"

    def __init__(self, name, conf, props, mesh):
        self.name = name
        myProps = props.getProps(name)
        myConf = conf.makeProps(name)

        self.type = myProps.get("type", "Constraints")
        self.conTable = myProps.get("conTable")

        myConf.set("type", self.type)
        myConf.set("conTable", self.conTable)

        if self.conTable is None:
            raise KeyError(self.__key__.format(self.type, "conTable"))

    def takeAction(self, action, globdat):
        if action == "GET_CONSTRAINTS":
            cons = globdat.get(self.conTable)
            loadScale = globdat.get("loadScale")
            raise NotImplementedError()
            # cons.scale?
        else:
            return False

#===========================================================================
#   LoadScaleModel
#===========================================================================


class LoadScaleModel(Model):
    """ The root of the model tree 

    Instance Members:
        name = model name
        type = model type ("LoadScale")
        model = child model

    Public Methods:
        MatrixModel(name, conf, props, mesh)
        takeAction(action, globdat)
    """

    def __init__(self, name, conf, props, mesh):
        """ Creates a node and its child """
        self.name = name
        myProps = props.getProps(name)
        myConf = conf.makeProps(name)

        self.type = myProps.get("type","LoadScale")
        myConf.set("type",self.type)

        # Create child
        self.model = ModelFactory("model", myConf, myProps, mesh)

    def takeAction(self, action, globdat):
        if action == "ADVANCE":
            self.model.takeAction(action, globdat)
            return True
