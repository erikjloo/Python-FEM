# Import Standard Libraries
from abc import ABCMeta, abstractmethod
from enum import IntEnum
import logging

#===========================================================================
#   Status
#===========================================================================


class Action(IntEnum):
    GET_MATRIX_0 = 0
    GET_EXT_VECTOR = 1
    GET_INT_VECTOR = 2
    GET_CONSTRAINTS = 3
    ADVANCE = 4
    COMMIT = 5


#===========================================================================
#   Model
#===========================================================================


class Model(metaclass=ABCMeta):
    """ Abstract Model Class
    
    Pure Virtual Methods:
        Model(name, conf, props, globdat)
        takeAction(action, globdat)
    
    Static Method:
        modelFactory(name, conf, props, globdat)
    """

    def __init__(self, name, conf, props, globdat):
        self.name = name

    def __del__(self):
        logging.debug("Cleaning %s model",self.name)

    @abstractmethod
    def takeAction(self, action, globdat):
        raise NotImplementedError()
    
    @staticmethod
    def modelFactory(name, conf, props, globdat):
        message = "Creating a {} model named {}"
        type = props.get("{}.type".format(name))
        logging.info(message.format(type, name))

        if type == "Matrix":
            return MatrixModel(name, conf, props, globdat)
        elif type == "Multi":
            return MultiModel(name, conf, props, globdat)
        elif type == "Solid":
            from solidModel import SolidModel
            return SolidModel(name, conf, props, globdat)
        elif type == "Truss":
            from trussModel import TrussModel
            return TrussModel(name, conf, props, globdat)
        elif type == "Periodic":
            from PBCmodel import PBCmodel
            return PBCmodel(name, conf, props, globdat)
        elif type == "PointLoad":
            return PointLoadModel(name, conf, props, globdat)
        elif type == "Constraints":
            return ConstraintsModel(name, conf, props, globdat)
        elif type == "LoadScale":
            return LoadScaleModel(name, conf, props, globdat)
        else:
            raise NotImplementedError("{} model not implemented".format(type))

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
        MatrixModel(name, conf, props, globdat)
        takeAction(action, globdat)
    """

    def __init__(self, name, conf, props, globdat):
        self.name = name
        myProps = props.getProps(name)
        myConf = conf.makeProps(name)

        self.type = myProps.get("type", "Matrix")
        myConf.set("type", self.type)

        # Create child
        self.model = self.modelFactory("model", myConf, myProps, globdat)

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
        MultiModel(name, conf, props, globdat)
        takeAction(action, globdat)
    """

    def __init__(self, name, conf, props, globdat):
        self.name = name
        myProps = props.getProps(name)
        myConf = conf.makeProps(name)

        self.type = myProps.get("type", "Multi")
        sub_models = myProps.get("models")

        myConf.set("type", self.type)
        myConf.set("models", sub_models)

        # Create children
        self.models = []
        for name in sub_models:
            model = self.modelFactory(name, myConf, myProps, globdat)
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
        PointLoadModel(name, conf, props, globdat)
        takeAction(action, globdat)
    """

    def __init__(self, name, conf, props, globdat):
        self.name = name
        myProps = props.getProps(name)
        myConf = conf.makeProps(name)

        self.type = myProps.get("type", "PointLoad")
        self.loadTable = myProps.get("loadTable")

        myConf.set("type", self.type)
        myConf.set("loadTable", self.loadTable)

        mesh = globdat.get("mesh")
        load = globdat.get(self.loadTable)
        self.rvals = load.initialize(mesh)

    def takeAction(self, action, globdat):
        if action == Action.GET_EXT_VECTOR:
            fext = globdat.get("fext")
            scale = globdat.get("loadScale")
            fext += scale*self.rvals
            globdat.set("fext",fext)
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
        ConstraintsModel(name, conf, props, globdat)
        takeAction(action, globdat)
    """

    def __init__(self, name, conf, props, globdat):
        self.name = name
        myProps = props.getProps(name)
        myConf = conf.makeProps(name)

        self.type = myProps.get("type", "Constraints")
        self.conTable = myProps.get("conTable")

        myConf.set("type", self.type)
        myConf.set("conTable", self.conTable)

        mesh = globdat.get("mesh")
        cons = globdat.get(self.conTable)
        self.rvals, self.sdof = cons.initialize(mesh)

    def takeAction(self, action, globdat):
        if action == Action.GET_CONSTRAINTS:
            cons = globdat.get(self.conTable)
            scale = globdat.find("loadScale",1)
            for idof in self.sdof:
                cons.addConstraint(idof, scale*self.rvals[idof])
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
        type = model type ("LoadScale")
        model = child model

    Public Methods:
        MatrixModel(name, conf, props, globdat)
        takeAction(action, globdat)
    """

    def __init__(self, name, conf, props, globdat):
        self.name = name
        myProps = props.getProps(name)
        myConf = conf.makeProps(name)

        self.type = myProps.get("type","LoadScale")
        myConf.set("type",self.type)

        # Create child
        self.model = self.modelFactory("model", myConf, myProps, globdat)

    def takeAction(self, action, globdat):
        if action == "ADVANCE":
            self.model.takeAction(action, globdat)
        elif action == Action.GET_CONSTRAINTS or action == Action.GET_EXT_VECTOR:
            # update thing,
            # call child
            return True
