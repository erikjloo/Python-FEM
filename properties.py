# Import Standard Libraries
import json
from pprint import pprint
from collections import defaultdict

#===========================================================================
#   Properties
#===========================================================================


class Properties(object):
    """ Properties

    Instance Members:
        properties = dictionary of properties
        
    Public Methods:
        Properties(my_dict=None)
        parseFile(file)
        writeFile(file)
        makeProps(props)
        dict = get(props)
        set(prop, value)
        props = getProps(props)
        print()
    """

    # Public:

    #-----------------------------------------------------------------------
    #   constructor
    #-----------------------------------------------------------------------

    def __init__(self, my_dict=None):
        """ Input: dictionary """
        if my_dict is None:
            self.properties = {} #defaultdict(dict)
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
        """ Input:  props = string of property name to be created
            Output: properties object (of shallow copy) of given props """
        self.properties[props] = {}
        return Properties(self.properties[props])

    #-----------------------------------------------------------------------
    #   getProps
    #-----------------------------------------------------------------------

    def getProps(self, props):
        """ Input: props = string of property names separated by '.'
            Output: properties object (of shallow copy) of given props """
        return Properties(self.get(props))

    #-----------------------------------------------------------------------
    #   get
    #-----------------------------------------------------------------------

    def get(self, props, default=None):
        """ Input:  props = string of property names separated by '.'
            Output: dictionary or value of properties of given props """
        props = props.split('.')

        try:
            tmp = self.properties[props[0]]
            for key in props[1:]:
                tmp = tmp[key]
            return tmp
        except KeyError:
            return default

    #-----------------------------------------------------------------------
    #   set
    #-----------------------------------------------------------------------

    def set(self, props, value):
        """ Input:  props = string of property names separated by '.'
                    value"""
        props = props.split('.')
                
        if len(props) is 1:
            self.properties[props[0]] = value
        elif len(props) is 2:
            if props[0] in self.properties:
                self.properties[props[0]][props[1]] = value
            else:
                self.properties[props[0]] = {}
                self.properties[props[0]][props[1]] = value
        else:
            print(" Cannot nest deeper than 2 ")

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
