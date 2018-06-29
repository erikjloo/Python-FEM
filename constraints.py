# Import Standard Libraries
import re
import scipy as np

#===========================================================================
#   Constraints
#===========================================================================


class Constraints(object):
    """ Constraints

    Static Members:
        __type_int__ = "Input idof is not int!"
        __type_int_list__ = "Input idofs is not int or list!"

    Instance Members:
        name = table name
        type = table type ("Constraints")
        ndof = number of degrees of freedom
        conspace = array of constraints per idof
        sdof = list of supported degrees of freedom
        
    Public Methods:
        Constraints(ndof=0)
        initialize(name, conf, props, mesh)
        readXML(path, mesh)
        addConstraint(idof, rval=0.0)
        addConstraints(idofs, rval=0.0)
        eraseConstraint(idof)
        eraseConstraints(idofs)
        disps = getDisps()
        ndof = dofCount()
        fdof = get_fdof()
        sdof = get_sdof()
    """

    # Static:
    __type_int__ = "Input inod is not int!"
    __type_int_list__ = "Input is not int or list!"

    # Public:

    #-----------------------------------------------------------------------
    #   constructor
    #-----------------------------------------------------------------------

    def __init__(self, ndof=0):
        """ Input: ndof = number of degrees of freedom """
        self.ndof = ndof
        self.sdof = []
        self.conspace = np.empty(ndof)
        self.conspace[:] = np.nan

    def resize(self, ndof):
        """ Input: ndof = new size external force vector """
        self.conspace.resize(ndof)

    #-----------------------------------------------------------------------
    #   initialize
    #-----------------------------------------------------------------------

    def initialize(self, name, conf, props, mesh):
        """ Input:  name = table name
                    conf = output properties
                    props = input properties
                    mesh = Mesh """
        self.name = name
        myProps = props.getProps(name)
        myConf = conf.makeProps(name)

        self.type = myProps.get("type", "Constraints")
        path = myProps.get("file")

        myConf.set("type", self.type)
        myConf.set("file", path)

        self.ndof = mesh.dofCount()
        self.conspace.resize(self.ndof)
        self.conspace[:] = np.nan
        self.readXML(path, mesh)
        print(path, "file read")

    #-----------------------------------------------------------------------
    #   readXML
    #-----------------------------------------------------------------------

    def readXML(self, path, mesh):
        """ Input:  path = path_to_file, mesh = Mesh """
        with open(path, 'r') as file:
            
            flag_c = False

            for line in file:
                if line.startswith("<Constraints>"):
                    flag_c = True
                elif line.startswith("</Constraints>"):
                    flag_c = False
                
                if flag_c is True and not line.startswith("<Constraints>"):
                    dof = re.findall(r"[a-zA-Z]+", line)[0]
                    [node, rval] = re.findall(r"[-+]?\d *\.\d+|[-+]?\d+", line)
                    print(" {}[{}] = {}".format(dof, node, rval))
                    idof = mesh.getDofIndex(int(node), dof)
                    self.addConstraint(idof, float(rval))

    #-----------------------------------------------------------------------
    #   Constraint Methods
    #-----------------------------------------------------------------------

    def addConstraint(self, idof, rval=0.0):
        """ Input:  idof = dof index
                    rval = prescribed displacement (default: 0.0) """
        if isinstance(idof, int):
            if idof not in self.sdof:
                self.conspace[idof] = rval
                self.sdof.append(idof)
        else:
            raise TypeError(self.__type_int__)
            
    def addConstraints(self, idofs, rval=0.0):
        """ Input:  idofs = (list of) dof indices
                    rval = prescribed displacement (default: 0.0) """
        if isinstance(idofs,(list,tuple,range,np.ndarray)):
            for idof in idofs:
                self.addConstraint(idof, rval)
        elif isinstance(idofs, int):
            self.addConstraint(idof, rval)
        else:
            raise TypeError(self.__type_int_list__)

    def eraseConstraint(self,  idof):
        """ Input: idof = prescribed dof index to be erased """
        if isinstance(idof, int):
            self.conspace[idof] = np.nan
            idx = self.sdof.index(idof)
            del self.sdof[idx]
        else:
            raise TypeError(self.__type_int__)

    def eraseConstraints(self, idofs):
        """ Input: idofs = (list of) prescribed dof indices to be erased """
        if isinstance(idofs, (list, tuple, range, np.ndarray)):
            for idof in idofs:
                self.eraseConstraint(idof)
        elif isinstance(idofs, int):
            self.eraseConstraint(idofs)
        else:
            raise TypeError(self.__type_int_list__)

    #-----------------------------------------------------------------------
    #   getDisps
    #-----------------------------------------------------------------------

    def getDisps(self):
        disps = np.zeros(self.ndof)
        disps[self.sdof] = self.conspace[self.sdof]
        return disps

    def getConstraints(self):
        return self.conspace[self.sdof]
        
    def get_fdof(self):
        fdof = np.argwhere(np.isnan(self.conspace)).transpose()[0]
        return fdof.tolist()

    def get_sdof(self):
        return self.sdof
