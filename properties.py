# Import Standard Libraries
import json
import scipy as np
from enum import IntEnum
from pprint import pprint


#===========================================================================
#   Properties
#===========================================================================


class Properties(object):

    __type_str__ = "Input is not str!"

    def __init__(self, my_dict=None):
        """ Input: dictionary """
        if my_dict is None:
            self.properties = {}
        else:
            self.properties = my_dict

    #-----------------------------------------------------------------------
    #   parseFile
    #-----------------------------------------------------------------------

    def parseFile(self, file):
        """ Input: file_path """
        with open(file, 'r') as f:
            self.properties = json.load(f)

    #-----------------------------------------------------------------------
    #   writeFile
    #-----------------------------------------------------------------------

    def writeFile(self, file):
        """ Input: file_path """
        with open(file, 'w') as f:
            json.dump(self.properties, f)

    #-----------------------------------------------------------------------
    #   makeProps
    #-----------------------------------------------------------------------

    def makeProps(self, props):
        pass

    #-----------------------------------------------------------------------
    #   get
    #-----------------------------------------------------------------------

    def get(self, props):
        """ Input: props = string of property names separated by '.'
            Output: dictionary of properties of given props """
        try:
            props = props.split('.')
        except:
            TypeError(self.__type_str__)

        if len(props) is 1:
            return self.properties[props[0]]
        if len(props) is 2:
            return self.properties[props[0]][props[1]]
        elif len(props) is 3:
            return self.properties[props[0]][props[1]][props[2]]
        else:
            print(" Cannot nest deeper than 3 ")

    #-----------------------------------------------------------------------
    #   set
    #-----------------------------------------------------------------------

    def set(self, prop, value):
        pass

    #-----------------------------------------------------------------------
    #   getProps
    #-----------------------------------------------------------------------

    def getProps(self, props):
        """ Input: props = string of property names separated by '.'
            Output: Properties object of given props """
        return Properties(self.get(props))

    def print(self):
        pprint(self.properties)


#===========================================================================
#   Example
#===========================================================================


if __name__ == "__main__":

    file = "Examples/2D_semicircle.pro"

    props = Properties()
    props.parseFile(file)
    props.print()

    matrix_props = props.getProps("model")
    matrix_props.print()
