# Import Standard Libraries
import re
import scipy as np
from warnings import warn

#===========================================================================
#   Constraints
#===========================================================================


class Constraints(object):
    """ Constraints

    Static Members:
        __type_int__ = "Input inod is not int!"
        __type_int_list__ = "Input is not int or list!"

    Instance Members:
        ndof = number of degrees of freedom
        conspace = array of constraints per idof
        sdof = list of supported degrees of freedom
        
    Public Methods:
        Constraints()
        initialize(props, mesh)
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
        self.ndof = ndof
        self.sdof = []
        self.conspace = np.empty(ndof)
        self.conspace[:] = np.nan

    #-----------------------------------------------------------------------
    #   initialize
    #-----------------------------------------------------------------------

    def initialize(self, props, mesh):
        """ Input:  props = Properties
                    mesh = Mesh """
        self.ndof = mesh.dofCount()
        self.conspace = np.empty(self.ndof)
        self.conspace[:] = np.nan

        try:
            myProps = props.getProps("input.constraints")
            path = myProps.get("file")
            self.readXML(path, mesh)
            print(path,"file read")
        except TypeError:
            warn(" No constraints provided ")

    #-----------------------------------------------------------------------
    #   readXML
    #-----------------------------------------------------------------------

    def readXML(self, path, mesh):
        """ Input:  path = path_to_file 
                    mesh = Mesh"""
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
        """ Input:  idofs = list of dof indices
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
            raise TypeError(self.__type_int_list__)

    def eraseConstraints(self, idofs):
        """ Input: idofs = list of prescribed dof indices to be erased """
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

    def dofCount(self):
        return self.ndof

    def get_fdof(self):
        fdof = np.argwhere(np.isnan(self.conspace)).transpose()[0]
        return fdof.tolist()

    def get_sdof(self):
        return self.sdof
