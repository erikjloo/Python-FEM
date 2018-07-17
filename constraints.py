# Import Standard Libraries
import re
import logging
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
        path = file path to constraints
        ndof = number of degrees of freedom
        rvals = array of constraints per idof
        sdof = list of supported degrees of freedom
        
    Public Methods:
        Constraints(name, conf, props)
        initialize(mesh)
        readXML(path, mesh)
        addConstraint(idof, rval=0.0)
        addConstraints(idofs, rval=0.0)
        eraseConstraint(idof)
        eraseConstraints(idofs)
        updateSolution(solu)
        ndof = dofCount()
        fdof = getFdof()
        sdof = getSdof()
    """

    # Static:
    __type_int__ = "Input inod is not int!"
    __type_int_list__ = "Input is not int or list!"

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

            self.type = myProps.get("type", "Constraints")
            self.path = myProps.get("file")

            myConf.set("type", self.type)
            myConf.set("file", self.path)

        elif isinstance(name, int):
            self.sdof = []
            self.rvals = np.empty(name)
            self.rvals[:] = np.nan

    #-----------------------------------------------------------------------
    #   initialize
    #-----------------------------------------------------------------------

    def initialize(self, mesh):
        self.sdof = []
        self.ndof = mesh.dofCount()
        self.rvals = np.empty(self.ndof)
        self.rvals[:] = np.nan
        self.readXML(self.path, mesh)
        logging.info("    %s file read",self.path)
        return self.rvals, self.sdof

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
                    [node, rval] = re.findall(r"[-+]?\d+ *\.\d+|[-+]?\d+", line)
                    logging.debug("    %s[%s] = %s",dof, node, rval)
                    idof = mesh.getDofIndex(int(node), dof)
                    self.addConstraint(idof, float(rval))

    #-----------------------------------------------------------------------
    #   Constraint Methods
    #-----------------------------------------------------------------------

    def addConstraint(self, idof, rval=0.0):
        """ Input:  idof = dof index
                    rval = prescribed displacement (default: 0.0) """
        if isinstance(idof, int):
            self.rvals[idof] = rval
            if idof not in self.sdof:
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
            self.rvals[idof] = np.nan
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
    #   updateSolution
    #-----------------------------------------------------------------------
    def updateSolution(self, solu):
        solu[self.sdof] = self.rvals[self.sdof]

    def getFdof(self):
        fdof = np.argwhere(np.isnan(self.rvals)).transpose()[0]
        return fdof.tolist()

    def getSdof(self):
        return self.sdof
