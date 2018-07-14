# Import Standard Libraries
import json
import inspect
from pprint import pprint

#===========================================================================
#   Properties
#===========================================================================


class Properties(object):
    """ Properties

    Static Members:
        __type_dict__ = " Key attribute is not a dict"

    Instance Members:
        properties = dictionary of properties
        
    Public Methods:
        Properties(my_dict=None)

        Parser/Writer:
            parseFile(file)
            writeFile(file)

        Prop Methods:
            *props = makeProps(props)
            *props = findProps(props)
            *props = getProps(props)
        
        Dict Methods:
            *dict/value = set(prop, value)
            *dict/value = find(props)
            *dict/value = get(props)
            *dict = getDict()

        Miscellaneous:
            print()
    """

    __type_dict__ = " {}'s value should be of type dict !"

    # Public:

    #-----------------------------------------------------------------------
    #   Constructor
    #-----------------------------------------------------------------------

    def __init__(self, my_dict=None):
        """ Input: dictionary """
        if my_dict is None:
            self.properties = {}
        else:
            self.properties = my_dict

    #-----------------------------------------------------------------------
    #   Parser/Writer
    #-----------------------------------------------------------------------

    def parseFile(self, file):
        """ Input: file_path """
        with open(file, 'r') as f:
            self.properties = json.load(f)

    def writeFile(self, file):
        """ Input:  file_path 
            Output: text file with configuration data """
        with open(file, 'w') as f:
            json.dump(self.properties, f, sort_keys=False, indent=4)

    #-----------------------------------------------------------------------
    #   Prop Methods
    #-----------------------------------------------------------------------

    def makeProps(self, props):
        """ Input:  props = string of property name to be created
            Output: properties object (of shallow copy) of given props """
        if props in self.properties:
            if isinstance(self.properties[props],dict):
                return Properties(self.properties[props])
        else:
            self.properties[props] = {}
            return Properties(self.properties[props])

    def findProps(self, props):
        """ Input: props = string of property name """
        if isinstance(self.find(props), dict):
            return Properties(self.find(props))
        else:
            raise TypeError(self.__type_dict__.format(props))

    def getProps(self, props):
        """ Input: props = string of property names separated by '.'
            Output: properties object (of shallow copy) of given props """
        if isinstance(self.get(props), dict):
            return Properties(self.get(props))
        else:
            raise TypeError(self.__type_dict__.format(props))

    #-----------------------------------------------------------------------
    #   Dict Methods
    #-----------------------------------------------------------------------

    def set(self, props, value):
        """ Input:  props = string of property names separated by '.'
            Output: value = value set to given props """
        props = props.split('.')
                
        if len(props) is 1:
            self.properties[props[0]] = value
            return self.properties[props[0]]
        elif len(props) is 2:
            if props[0] in self.properties:
                self.properties[props[0]][props[1]] = value
            else:
                self.properties[props[0]] = {}
                self.properties[props[0]][props[1]] = value
            return self.properties[props[0]][props[1]]
        else:
            raise ValueError(" Cannot nest deeper than 2 ")

    def find(self, props, default=None):
        """ Input:  props = string of property names separated by '.'
            Output: value = values of given props or default """
        props = props.split('.')

        try:
            tmp = self.properties[props[0]]
            for key in props[1:]:
                tmp = tmp[key]
            return tmp
        except KeyError:
            return default

    def get(self, props, default=None):
        """ Input:  props = string of property names separated by '.'
            Output: value = values of given props or KeyError if None """
        tmp = self.find(props, default)

        if tmp is None:
            stack = inspect.stack()
            name = stack[1][0].f_locals["self"].__class__
            raise KeyError("{}: {} not specified !".format(name, props))
        else:
            return tmp

    def getDict(self):
        """ Output: dictionary of properties """
        return self.properties

    #-----------------------------------------------------------------------
    #   Miscellaneous
    #-----------------------------------------------------------------------

    def print(self):
        pprint(self.properties)
