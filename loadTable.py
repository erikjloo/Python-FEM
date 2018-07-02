# Import Standard Libraries
import re
import scipy as np

#===========================================================================
#   LoadTable
#===========================================================================


class LoadTable(object):
    """ LoadTable

    Static Members:
        __type_int__ = "Input idof is not int!"
        __type_int_list__ = "Input idofs is not int or list!"

    Instance Members:
        name = table name
        type = table type ("Loads")
        path = file path to loads
        rvals = array of point loads per idof

    Public Methods:
        LoadTable(name, conf, props)
        initialize(mesh)
        readXML(path, mesh)
        setLoad(idof, rval)
        setLoads(idofs, rval)
        addLoad(idof, rval)
        addLoads(idofs, rval)
        rvals = getLoads()
    """

    # Static:
    __type_int__ = "Input idof is not int!"
    __type_int_list__ = "Input idofs is not int or list!"

    # Public:

    #-----------------------------------------------------------------------
    #   constructor
    #-----------------------------------------------------------------------

    def __init__(self, name, conf=None, props=None):
        """ Input:  name = table name or ndof
                    conf = output properties
                    props = input properties """
        if isinstance(name, str):
            self.name = name
            myProps = props.getProps(name)
            myConf = conf.makeProps(name)

            self.type = myProps.get("type", "Loads")
            self.path = myProps.get("file")

            myConf.set("type", self.type)
            myConf.set("file", self.path)

        elif isinstance(name, int):
            self.rvals = np.empty(name)
            self.rvals[:] = 0

    #-----------------------------------------------------------------------
    #   initialize
    #-----------------------------------------------------------------------

    def initialize(self, mesh):
        self.rvals = np.zeros(mesh.dofCount())
        self.readXML(self.path, mesh)
        print(self.path, "file read")
        return self.rvals

    #-----------------------------------------------------------------------
    #   readXML
    #-----------------------------------------------------------------------

    def readXML(self, path, mesh):
        """ Input:  path = path_to_file, mesh = Mesh """
        with open(path, 'r') as file:

            flag_c = False

            for line in file:
                if line.startswith("<Loads>"):
                    flag_c = True
                elif line.startswith("</Loads>"):
                    flag_c = False

                if flag_c is True and not line.startswith("<Loads>"):
                    dof = re.findall(r"[a-zA-Z]+", line)[0]
                    [node, rval] = re.findall(r"[-+]?\d *\.\d+|[-+]?\d+", line)
                    print(" {}[{}] = {}".format(dof, node, rval))
                    idof = mesh.getDofIndex(int(node), dof)
                    self.addLoad(idof, float(rval))

    #-----------------------------------------------------------------------
    #   Load Methods
    #-----------------------------------------------------------------------
    
    def setLoad(self,  idof, rval):
        """ Input: idof = dof index, rval = load to be set """
        if isinstance(idof, int):
            self.rvals[idof] = rval
        else:
            raise TypeError(self.__type_int__)

    def setLoads(self, idofs, rval):
        """ Input: idofs = (list of) dof indices, rval = load to be set"""
        if isinstance(idofs, (int,list)):
            self.rvals[idofs] = rval
        else:
            raise TypeError(self.__type_int_list__)

    def addLoad(self, idof, rval):
        """ Input: idof = dof index, rval = load to be added """
        if isinstance(idof, int):
            self.rvals[idof] += rval
        else:
            raise TypeError(self.__type_int__)

    def addLoads(self, idofs, rval):
        """ Input: idofs = (list of) dof indices, rval = load to be added """
        if isinstance(idofs, (int,list)):
            self.rvals[idofs] += rval
        else:
            raise TypeError(self.__type_int_list__)

    #-----------------------------------------------------------------------
    #   getLoads
    #-----------------------------------------------------------------------

    def getLoads(self):
        return self.rvals
