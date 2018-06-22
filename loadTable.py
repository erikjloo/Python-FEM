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
        __type_int__ = "Input inod is not int!"
        __type_int_list__ = "Input is not int or list!"

    Instance Members:
        ndof = number of degrees of freedom
        loads = array of point loads per idof

    Public Methods:
        LoadTable(ndof=0)
        initialize(props, mesh)
        readXML(path, mesh)
        setLoad(idof, rval)
        setLoads(idofs, rval)
        addLoad(idof, rval)
        addLoads(idofs, rval)
        loads = getLoads()
    """

    # Static:
    __type_int__ = "Input inod is not int!"
    __type_int_list__ = "Input is not int or list!"

    # Public:

    #-----------------------------------------------------------------------
    #   constructor
    #-----------------------------------------------------------------------

    def __init__(self, ndof=0):
        self.ndof = ndof
        self.loads = np.empty(ndof)
        self.loads[:] = 0

    #-----------------------------------------------------------------------
    #   initialize
    #-----------------------------------------------------------------------

    def initialize(self, props, mesh):
        """ Input:  props = properties """

        self.ndof = mesh.dofCount()
        self.loads = np.empty(self.ndof)
        self.loads[:] = 0

        try:
            myProps = props.getProps("input.loads")
            path = myProps.get("file")
            self.readXML(path, mesh)
            print(path, "file read")
        except TypeError:
            warn(" No loads provided ")

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
        if isinstance(idof, int):
            self.loads[idof] = rval
        else:
            raise TypeError(self.__type_int__)

    def setLoads(self, idofs, rval):
        if isinstance(idofs, list):
            self.loads[idofs] = rval
        elif isinstance(idofs, int):
            self.setLoad(idofs,rval)
        else:
            raise TypeError(self.__type_int_list__)

    def addLoad(self, idof, rval):
        if isinstance(idof, int):
            self.loads[idof] += rval
        else:
            raise TypeError(self.__type_int__)

    def addLoads(self, idofs, rval):
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
