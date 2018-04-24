# Import Local Libraries
from properties import Properties
from solidModel import SolidModel
from PBCmodel import PBCmodel


#===========================================================================
#   Example
#===========================================================================


class ModelFactory(object):

    def __init__(self, props, name="model"):
        self.type = props.get(name+".type")
        self.props = props
        self.name = name

    def createModel(self):

        if self.type == "Multi":
            models = MultiModel(self.props)
            return models.createModels()

        elif self.type == "Solid":

            return SolidModel(self.props, self.name)

        elif self.type == "Periodic":

            return PBCmodel(self.props, self.name)


#===========================================================================
#   MultiModel
#===========================================================================


class MultiModel(ModelFactory):

    def __init__(self, props):
        self.names = props.get("model.models")
        self.props = props.getProps("model")

    def createModels(self):

        models = []

        for name in self.names:
            
            sub_model = ModelFactory(self.props, name)
            sub_model = sub_model.createModel()
            models.append(sub_model) 
        
        return models


#===========================================================================
#   Example
#===========================================================================


if __name__ == "__main__":

    file = "Examples/square.pro"

    props = Properties()
    props.parseFile(file)

    model = ModelFactory(props)
    models = model.createModel()
    
    models[0].initialize()
    ndof = models[0].dofCount()
    print(ndof)
