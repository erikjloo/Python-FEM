from properties import Properties

class MultiModel(object):

    def __init__(self, props):
        props = props.getProps("model")
        self.models = props.get("models")
        self.createModels(props)

    def createModels(self, props):

        for model in enumerate(self.models):
            print(model+".type")
            type = props.getProps(model+".type")
            model = ModelFactory(props)

class ModelFactory(object):
    
    def __init__(self, props):
        self.type = props.get("model.type")

    def createModel(self, type, props, globdat=None):
        if self.type == "Multi":
            return MultiModel(props) 
        elif self.type == "Solid":
            from solidModel import SolidModel
            return SolidModel("model", props)
        elif self.type == "Periodic":
            from PBCmodel import PBCmodel
            return PBCmodel("model", props)


#===========================================================================
#   Example
#===========================================================================


if __name__ == "__main__":

    file = "Examples/square.pro"

    props = Properties()
    props.parseFile(file)

    model = ModelFactory()
    model.createModel(props)
