# Import Standard Libraries
import scipy as np
import warnings


#===========================================================================
#   DofSpace
#===========================================================================


class DofSpace(object):
    """ Dof Space 
    
    Static Members:
        __type__ = "Input is not list or array!"
        __type_int__ = "Input is not int!"
        __type_str__ = "Input is not str!"
        __renumber__ = "Erasing dofs: Dof numbers will be renumbered!"

    Instance Members:
        nrow = number of rows (nodes)
        types = list of dof type names
        dofspace = array of dof indices (idofs)
            i row corresponds to inod
            j column corresponds to jtype
            idof = dofspace[inod,jtype]
        ndof = last dof index

    Public Methods:
        DofSpace(nrow, ntyp)

        Row Methods:
            nrow = addRow()
            nrow = addRows(nrow)
            eraseRow(irow)
            eraseRows(irows)
            nrow = rowCount()
            
        Type Methods:
            addType(dof)
            addTypes(dofs)
            setType(jtype, dof)
            eraseType(dof)
            eraseTypes(dofs)
            ntyp = typeCount()
            dof = getTypeName(jtype)

        Dof Methods:
            addDof(inod, dofs)
            addDofs(inodes, dofs)
            eraseDof(inod, dofs)
            eraseDofs(inodes, dofs)
            ndof = dofCount()
            idof = getDofIndex(inod, dof)
            idofs = getDofIndices(inodes, dofs)
            printDofSpace()

    Private Methods:
        __renumberDofs()
    """
    
    # Static:
    __type__ = "Input is not list or array!"
    __type_int__ = "Input is not int!"
    __type_str__ = "Input is not str!"
    __renumber__ = "Erasing dofs: Dof numbers will be renumbered!"

    # Public:
    def __init__(self, nrow, ntyp=1):
        """ Input: nrow = number of rows (nodes), ntyp = number of dof types """
        self.nrow = nrow
        self.types = []
        self.dofspace = np.empty((nrow, ntyp))
        self.dofspace[:] = np.nan
        self.ndof = 0

    #-------------------------------------------------------------------
    #   Row Methods
    #-------------------------------------------------------------------

    def addRow(self):
        """ Adds a new row to dofspace """
        r_reqd = np.size(self.dofspace, 1)
        r_new = np.empty((1,r_reqd))
        r_new[:] = np.nan
        self.dofspace = np.vstack((self.dofspace, r_new))
        self.nrow += 1
        return self.nrow

    def addRows(self, nrow):
        """ Input: nrow = number of rows (nodes) to be added """
        r_reqd = np.size(self.dofspace, 1)
        r_new = np.empty((nrow,r_reqd))
        r_new[:] = np.nan
        self.dofspace = np.vstack((self.dofspace, r_new))
        self.nrow += nrow
        return self.nrow
    
    def eraseRow(self, irow):
        """ Input: irow = index of row (node) to be erased """
        self.dofspace = np.delete(self.dofspace, irow, 0)
        warnings.warn(self.__renumber__)
        self.__renumberDofs()
        self.nrow -= 1
    
    def eraseRows(self, irows):
        """ Input: irows = indices of rows (nodes) to be erased """
        if isinstance(irows, (list, tuple, range, np.ndarray)):
            for irow in sorted(irows, reverse=True):
                self.eraseRow(irow)
        else:
            self.eraseRow(irows)

    def rowCount(self):
        """ Output: number of rows (nodes) """
        return self.nrow

    #-------------------------------------------------------------------
    #   Type Methods
    #-------------------------------------------------------------------

    def addType(self, dof):
        """ Input: dof = string of dof type name """
        if isinstance(dof, str):
            self.types.append(dof)
        else:
            raise TypeError(self.__type_str__)
        # Check if dofspace has enough columns for all dof types
        c_reqd = len(self.types) - np.size(self.dofspace, 1)
        if c_reqd > 0:
            c_new = np.empty((self.nrow, c_reqd))
            c_new[:] = np.nan
            self.dofspace = np.hstack((self.dofspace, c_new))

    def addTypes(self, dofs):
        """ Input: dofs =  list of string of dof type names """
        if isinstance(dofs, (list,np.ndarray)):
            for dof in dofs:
                self.addType(dof)
        else:
            raise TypeError(self.__type__)

    def setType(self, jtype, dof):
        """ Input: jtype = dof type index, dof = string of dof type name """
        if isinstance(dof, str):
            self.types[jtype] = dof
        else:
            raise TypeError(self.__type_str__)

    def eraseType(self, dof):
        """ Input: dof = string of dof type name to be erased """
        # Delete column from dofspace
        jtype = self.types.index(dof)
        self.dofspace = np.delete(self.dofspace, jtype, 1)
        # Delete type
        del self.types[jtype]
        # Renumber remaining dofs
        warnings.warn(self.__renumber__)
        self.__renumberDofs()

    def eraseTypes(self, dofs):
        """ Input: dofs = list of strings of dof type names to be erased """
        for dof in dofs:
            # Delete column from dofspace
            jtype = self.types.index(dof)
            self.dofspace = np.delete(self.dofspace, jtype, 1)
            # Delete type
            del self.types[jtype]
        # Renumber remaining dofs
        warnings.warn(self.__renumber__)
        self.__renumberDofs()

    def typeCount(self):
        """ Output: number of dof types """
        return len(self.types)

    def getTypeName(self, jtype):
        """ Input: jtype = dof type index
            Output: string of dof type name """
        return self.types[jtype]

    #-------------------------------------------------------------------
    #   Dof Methods
    #-------------------------------------------------------------------

    def addDof(self, inod, dofs):
        """ Input: inod = node index, dof = (list of) strings of dof names """
        # Need to check if inod is int or
        # 2 or > dofs will have same idof
        if isinstance(inod, int):
            for dof in dofs:
                jtype = self.types.index(dof)
                if np.isnan(self.dofspace[inod, jtype]):
                    self.dofspace[inod, jtype] = self.ndof
                    self.ndof += 1
        else:
            raise TypeError(self.__type_int__)

    def addDofs(self, inodes, dofs):
        """ Input: inodes = node indices, dofs = (list of) strings of dof names """
        if isinstance(inodes, (list,tuple,range,np.ndarray)):
            for inod in inodes:
                self.addDof(inod, dofs)
        else:  # addDof checks if inodes is an int
            self.addDof(inodes, dofs)

    def eraseDof(self, inod, dofs):
        """ Input: inod = node index, dofs = (list of) strings of dof names """
        for dof in dofs:
            jtype = self.types.index(dof)
            self.dofspace[inod, jtype] = np.nan
        warnings.warn(self.__renumber__)
        self.__renumberDofs()

    def eraseDofs(self, inodes, dofs):
        """ Input: inodes = node indices, dofs = (list of) strings of dof names """
        self.eraseDof(inodes, dofs)

    def dofCount(self):
        """ Output: ndof = number of degrees of freedom """
        return self.ndof

    def getDofIndex(self, inod, dof):
        """ Input: inod = node index, dof = string of dof name
            Output: idof = dof index """
        jtype = self.types.index(dof)
        return int(self.dofspace[inod,jtype])

    def getDofIndices(self, inodes, dofs):
        """ Input: inodes = node indices, dofs = list of strings of dof names
            Output: idofs = list of dof indices """
        idofs = []
        if isinstance(inodes, (list,tuple,range,np.ndarray)):
            for inod in inodes:
                for dof in dofs:
                    idofs.append(self.getDofIndex(inod, dof))
        else:
            for dof in dofs:
                idofs.append(self.getDofIndex(inodes, dof))
        return idofs

    def printDofSpace(self):
        """ Prints the dof space with node numbers and dof type names """
        print("\n", self.types)
        i = 0
        for row in self.dofspace:
            print(i, row)
            i += 1
        print("\n")

    # Private:
    def __renumberDofs(self):
        """ Renumbers all defined dofs from 0 to ndof """
        self.ndof = 0
        for i in range(np.size(self.dofspace, 0)):
            for j in range(np.size(self.dofspace, 1)):
                if ~np.isnan(self.dofspace[i, j]):
                    self.dofspace[i, j] = self.ndof
                    self.ndof += 1


#===========================================================================
#   Example
#===========================================================================


if __name__ == '__main__':

    print("\nDofSpace :")
    dofs = DofSpace(1, 1)
    dofs.addRows(10)
    dofs.eraseRows([5,6,7,8,9,10])
    dofs.addRow()

    dofs.addType('u')
    dofs.addTypes(['v', 'j', 'rotx', 'roty', 'rotz'])
    dofs.setType(2, 'w')
    dofs.addDofs((0,2,4,5), ['u', 'v', 'w', 'rotx'])
    dofs.printDofSpace()

    dofs.eraseType('w')
    dofs.printDofSpace()
    dofs.eraseTypes(['rotx', 'roty', 'rotz'])
    dofs.printDofSpace()

    dofs.addType('w')
    dofs.addDof(0,'u')
    dofs.addDof(0,['u','v','w'])
    dofs.addDofs(range(6),['u','v','w'])
    dofs.printDofSpace()

    dofs.eraseDof(0, ['u', 'v','w'])
    dofs.eraseDofs([1,2],'u')
    dofs.eraseDofs(1,'v')

    dofs.eraseRow(4)
    dofs.printDofSpace()
    print("Dof count :", dofs.dofCount())
