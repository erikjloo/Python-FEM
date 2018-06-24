# Import Standard Libraries
import re
import scipy as np
from warnings import warn

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
        loads = array of point loads per idof

    Public Methods:
        LoadTable(ndof=0)
        initialize(name, conf, props, mesh)
        readXML(path, mesh)
        setLoad(idof, rval)
        setLoads(idofs, rval)
        addLoad(idof, rval)
        addLoads(idofs, rval)
        loads = getLoads()
    """

    # Static:
    __type_int__ = "Input idof is not int!"
    __type_int_list__ = "Input idofs is not int or list!"

    # Public:

    #-----------------------------------------------------------------------
    #   constructor
    #-----------------------------------------------------------------------

    def __init__(self, ndof=0):
        """ Input: ndof = number of degrees of freedom """
        self.loads = np.empty(ndof)
        self.loads[:] = 0

    def resize(self, ndof):
        """ Input: ndof = new size external force vector """
        self.loads.resize(ndof)

    #-----------------------------------------------------------------------
    #   initialize
    #-----------------------------------------------------------------------

    def initialize(self, name, conf, props, mesh):
        """ Input:  name = table name
                    conf = output properties
                    props = input properties
                    mesh = Mesh """
        self.resize(mesh.dofCount())
        self.loads[:] = 0
        self.name = name
        myProps = props.getProps(name)
        myConf = conf.makeProps(name)

        self.type = myProps.get("type","Loads")
        path = myProps.get("file")

        myConf.set("type", self.type)
        myConf.set("file",path)

        self.readXML(path, mesh)
        print(path, "file read")

    #-----------------------------------------------------------------------
    #   readXML
    #-----------------------------------------------------------------------

    def readXML(self, path, mesh):
        """ Input: path = path_to_file """
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
            self.loads[idof] = rval
        else:
            raise TypeError(self.__type_int__)

    def setLoads(self, idofs, rval):
        """ Input: idofs = (list of) dof indices, rval = load to be set"""
        if isinstance(idofs, list):
            self.loads[idofs] = rval
        elif isinstance(idofs, int):
            self.setLoad(idofs,rval)
        else:
            raise TypeError(self.__type_int_list__)

    def addLoad(self, idof, rval):
        """ Input: idof = dof index, rval = load to be added """
        if isinstance(idof, int):
            self.loads[idof] += rval
        else:
            raise TypeError(self.__type_int__)

    def addLoads(self, idofs, rval):
        """ Input: idofs = (list of) dof indices, rval = load to be added """
        if isinstance(idofs, list):
            self.loads[idofs] += rval
        elif isinstance(idofs, int):
            self.addLoad(idofs, rval)
        else:
            raise TypeError(self.__type_int_list__)

    #-----------------------------------------------------------------------
    #   getLoads
    #-----------------------------------------------------------------------

    def getLoads(self):
        return self.loads
