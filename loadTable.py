# Import Standard Libraries
import re
import scipy as np

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
        LoadTable()
        initialize(props, mesh)
        readXML(path, mesh)
        addLoad(idof, rval=0.0)
        addLoads(idofs, rval=0.0)
        loads = getLoads()
    """

    # Static:
    __type_int__ = "Input inod is not int!"
    __type_int_list__ = "Input is not int or list!"

    # Public:

    #-----------------------------------------------------------------------
    #   constructor
    #-----------------------------------------------------------------------

    def __init__(self):
        self.ndof = 0
        self.loads = np.empty(0)
        self.loads[:] = np.nan

    #-----------------------------------------------------------------------
    #   initialize
    #-----------------------------------------------------------------------

    def initialize(self, props, mesh):
        """ Input:  props = properties """

        self.ndof = mesh.dofCount()
        self.loads = np.empty(self.ndof)
        self.loads[:] = 0

        try:
            props = props.getProps("input.loads")
            if props.get("type") == "Input":
                path = props.get("file")
                self.readXML(path, mesh)
                print(path, "file read")
        except:
            KeyError(" No Loads provided ")

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
        """ Input: idof = prescribed dof index to be erased """
        if isinstance(idof, int):
            self.loads[idof] = rval
        else:
            raise TypeError(self.__type_int_list__)

    def setLoads(self, idofs, rval):
        """ Input: idofs = list of prescribed dof indices to be erased """
        if isinstance(idofs, list):
            self.loads[idofs] = rval
        else: # setLoad checks if idof is int
            self.setLoad(idofs,rval)

    def addLoad(self, idof, rval):
        if isinstance(idof, int):
            self.loads[idof] += rval
        else:
            raise TypeError(self.__type_int_list__)

    def addLoads(self, idofs, rval):
        if isinstance(idofs, list):
            self.loads[idofs] += rval
        else: # addLoad checks if idofs is int
            self.addLoad(idofs, rval)

    #-----------------------------------------------------------------------
    #   getLoads
    #-----------------------------------------------------------------------

    def getLoads(self):
        return self.loads
