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
        loads = array of point loads per idof

    Public Methods:
        LoadTable(name, conf, props, mesh)
        initialize(mesh)
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

    def __init__(self, name, conf=None, props=None, mesh=None):
        """ Input:  name = table name or ndof
                    conf = output properties
                    props = input properties
                    mesh = Mesh """
        if isinstance(name, str):
            self.name = name
            myProps = props.getProps(name)
            myConf = conf.makeProps(name)

            self.type = myProps.get("type", "Loads")
            self.path = myProps.get("file")

            myConf.set("type", self.type)
            myConf.set("file", self.path)

        elif isinstance(name, int):
            self.loads = np.empty(name)
            self.loads[:] = 0

    #-----------------------------------------------------------------------
    #   initialize
    #-----------------------------------------------------------------------

    def initialize(self, mesh):
        self.loads = np.zeros(mesh.dofCount())
        self.readXML(self.path, mesh)
        print(self.path, "file read")

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
