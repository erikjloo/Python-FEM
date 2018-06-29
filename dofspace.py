# Import Standard Libraries
import scipy as np
from warnings import warn

#===========================================================================
#   DofSpace
#===========================================================================


class DofSpace(object):
    """ Dof Space 
    
    Static Members:
        __type_int__ = "Input is not int!"
        __type_str__ = "Input is not str!"
        __type_int_list__ = "Input is not int or list!"
        __type_str_list__ = "Input is not list or str!"
        __type_dof__ = "Input inod is not int or dof is not str!"
        __renumber__ = "Erasing dofs: Dof numbers will be renumbered!"
        
    Instance Members:
        nrow = number of rows (nodes)
        types = list of dof names
        dofspace = array of dof indices (idofs)
            i row corresponds to inod
            j column corresponds to jtype
            idof = dofspace[inod,jtype]
        ndof = number of degrees of freedom

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
            idof[s] = getDofIndex(inod, dof)
            idofs = getDofIndices(inodes, dofs)
        
        Miscellaneous:
            printDofSpace(irows)

    Private Methods:
        __renumberDofs()
    """
    
    # Static:
    __type_int__ = "Input is not int!"
    __type_str__ = "Input is not str!"
    __type_int_list__ = "Input is not int or list!"
    __type_str_list__ = "Input is not str or list!"
    __type_dof__ = "Input inod is not int or dof is not str!"
    __renumber__ = "Erasing dofs: Dof numbers will be renumbered!"

    # Public:

    #-------------------------------------------------------------------
    #   constructor
    #-------------------------------------------------------------------

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
        ntyp = np.size(self.dofspace, 1)
        r_new = np.empty((1, ntyp))
        r_new[:] = np.nan
        self.dofspace = np.vstack((self.dofspace, r_new))
        self.nrow += 1
        return self.nrow

    def addRows(self, nrow):
        """ Input: nrow = number of rows (nodes) to be added """
        ntyp = np.size(self.dofspace, 1)
        r_new = np.empty((nrow, ntyp))
        r_new[:] = np.nan
        self.dofspace = np.vstack((self.dofspace, r_new))
        self.nrow += nrow
        return self.nrow
    
    def eraseRow(self, irow):
        """ Input: irow = index of row (node) to be erased """
        if isinstance(irow, int):
            self.dofspace = np.delete(self.dofspace, irow, 0)
            warn(self.__renumber__)
            self.__renumberDofs()
            self.nrow -= 1
        else:
            raise TypeError(self.__type_int__)
    
    def eraseRows(self, irows):
        """ Input: irows = (list of) indices of rows (nodes) to be erased """
        if isinstance(irows, (list, tuple, range, np.ndarray)):
            for irow in sorted(irows, reverse=True):
                self.eraseRow(irow)
        elif isinstance(irows, int):
            self.eraseRow(irows)
        else:
            raise TypeError(self.__type_int_list__)

    def rowCount(self):
        """ Output: number of rows (nodes) """
        return self.nrow

    #-------------------------------------------------------------------
    #   Type Methods
    #-------------------------------------------------------------------

    def addType(self, dof):
        """ Input: dof = string of dof name """
        if isinstance(dof, str):
            if dof not in self.types:
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
        """ Input: dofs =  (list of) strings of dof names """
        if isinstance(dofs, (list, tuple, np.ndarray)):
            for dof in dofs:
                self.addType(dof)
        elif isinstance(dofs, str):
            self.addType(dofs)
        else:
            raise TypeError(self.__type_str_list__)

    def setType(self, jtype, dof):
        """ Input: jtype = dof type index, dof = string of dof name """
        if isinstance(dof, str):
            self.types[jtype] = dof
        else:
            raise TypeError(self.__type_str__)

    def eraseType(self, dof):
        """ Input: dof = string of dof name to be erased """
        if isinstance(dof, str):
            # Delete column from dofspace
            jtype = self.types.index(dof)
            self.dofspace = np.delete(self.dofspace, jtype, 1)
            # Delete type
            del self.types[jtype]
            # Renumber remaining dofs
            warn(self.__renumber__)
            self.__renumberDofs()
        else:
            raise TypeError(self.__type_str__)

    def eraseTypes(self, dofs):
        """ Input: dofs = (list of) strings of dof names to be erased """
        if isinstance(dofs, (list, tuple, np.ndarray)):
            for dof in dofs:
                # Delete column from dofspace
                jtype = self.types.index(dof)
                self.dofspace = np.delete(self.dofspace, jtype, 1)
                # Delete type
                del self.types[jtype]
            # Renumber remaining dofs
            warn(self.__renumber__)
            self.__renumberDofs()
        elif isinstance(dofs, str):
            self.eraseType(dofs)
        else:
            raise TypeError(self.__type_str_list__)

    def typeCount(self):
        """ Output: number of dof types """
        return len(self.types)

    def getTypeName(self, jtype):
        """ Input: jtype = dof type index
            Output: string of dof name """
        return self.types[jtype]

    #-------------------------------------------------------------------
    #   Dof Methods
    #-------------------------------------------------------------------

    def addDof(self, inod, dofs):
        """ Input: inod = node index, dofs = (list of) strings of dof names """
        # Need to check if inod is int or
        # 2 or > dofs will have same idof
        if isinstance(dofs, (list, tuple, np.ndarray)):
            for dof in dofs:
                self.addDof(inod,dof)
        elif isinstance(dofs, str) and isinstance(inod, int):
            jtype = self.types.index(dofs)
            if np.isnan(self.dofspace[inod, jtype]):
                self.dofspace[inod, jtype] = self.ndof
                self.ndof += 1
        else:
            raise TypeError(self.__type_dof__)

    def addDofs(self, inodes, dofs):
        """ Input: inodes = (list of) node indices, dofs = (list of) strings of dof names """
        if isinstance(inodes, (list, tuple, range, np.ndarray)):
            for inod in inodes:
                self.addDof(inod, dofs)
        elif isinstance(inodes, int):
            self.addDof(inodes, dofs)
        else:
            raise TypeError(self.__type_int_list__)

    def eraseDof(self, inod, dofs):
        """ Input: inodes = (list of) node indices, dofs = (list of) strings of dof names """
        if isinstance(dofs, (list, tuple, np.ndarray)):
            for dof in dofs:
                jtype = self.types.index(dof)
                self.dofspace[inod, jtype] = np.nan
        elif isinstance(dofs, str):
            jtype = self.types.index(dofs)
            self.dofspace[inod, jtype] = np.nan
        else:
            raise TypeError(self.__type_str_list__)
        warn(self.__renumber__)
        self.__renumberDofs()

    def eraseDofs(self, inodes, dofs):
        """ Input: inodes = (list of) node indices, dofs = (list of) strings of dof names """
        self.eraseDof(inodes, dofs)

    def dofCount(self):
        """ Output: ndof = number of degrees of freedom """
        return self.ndof

    def getDofIndex(self, inod, dofs=None):
        """ Input: inod = node index, dofs = (list of) strings of dof names
            Output: idof = (list of) dof indices """
        idofs = []
        if isinstance(dofs, (list, tuple, np.ndarray)):
            for dof in dofs:
                idofs.append(self.getDofIndex(inod, dof))
        elif isinstance(dofs, str) and isinstance(inod, int):
            jtype = self.types.index(dofs)
            idofs = int(self.dofspace[inod, jtype])
        elif dofs is None and isinstance(inod, int):
            for dof in self.dofspace[inod]:
                if not np.isnan(dof):
                    idofs.append(int(dof))
        else:
            raise TypeError(self.__type_dof__)

        return idofs

    def getDofIndices(self, inodes, dofs=None):
        """ Input: inodes = (list of) node indices, dofs = (list of) strings of dof names
            Output: idofs = (list of) dof indices """
        if isinstance(inodes, (list, tuple, range, np.ndarray)):
            idofs = []
            for inod in inodes:
                jdofs = self.getDofIndex(inod, dofs)
                if isinstance(jdofs, int):
                    idofs.append(jdofs)
                else:
                    idofs += jdofs
        elif isinstance(inodes, int):
            idofs = self.getDofIndex(inodes, dofs)
        else:
            raise TypeError(self.__type_int_list__)
        return idofs

    #-------------------------------------------------------------------
    #   Miscellaneous
    #-------------------------------------------------------------------

    def printDofSpace(self, rows=None):
        """ Input: rows = (list of) row indices to be printed """
        """ Prints dof space with (given) node numbers and dof names """
        print("\n", self.types)
        if rows is None:
            for i, row in enumerate(self.dofspace):
                print(i, row)
            print("\n")
        elif isinstance(rows, (list, tuple, range, np.ndarray)):
            for i, row in zip(rows, self.dofspace[rows,:]):
                print(i, row)
        elif isinstance(rows, int):
            print(rows, self.dofspace[rows,:])
        else:
            raise TypeError(self.__type_int_list__)

    # Private:
    def __renumberDofs(self):
        """ Renumbers all defined dofs from 0 to ndof-1 """
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
    dofs.eraseRows([5, 6, 7, 8, 9, 10])
    dofs.addRow()

    dofs.addType('u')
    dofs.addTypes(['v', 'j', 'rotx', 'roty', 'rotz'])
    dofs.setType(2, 'w')
    dofs.addDofs((0, 2, 4, 5), ['u', 'v', 'w', 'rotx'])
    dofs.printDofSpace()

    dofs.eraseType('w')
    dofs.printDofSpace()
    dofs.eraseTypes(['rotx', 'roty', 'rotz'])
    dofs.printDofSpace()

    dofs.addType('w')
    dofs.addDofs([0,1], 'u')
    dofs.addDof(0, ['u', 'v', 'w'])
    dofs.addDofs(range(6), ['u', 'v', 'w'])
    dofs.printDofSpace()

    dofs.eraseDof(0, ['u', 'v', 'w'])
    dofs.eraseDofs([1, 2], 'u')
    dofs.eraseDofs(1, 'v')

    dofs.eraseRow(4)
    dofs.printDofSpace()
    print("Dof count :", dofs.dofCount())

    idofs = dofs.getDofIndices(3)
    idofs = dofs.getDofIndices(3, "u")
    idofs = dofs.getDofIndices(3, ["u", "v"])

    idofs = dofs.getDofIndices([3, 4])
    idofs = dofs.getDofIndices([3, 4], "u")
    idofs = dofs.getDofIndices([3, 4], ["u","v"])
    print("idofs =", idofs)
