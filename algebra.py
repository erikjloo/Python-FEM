# Import Standard Libraries
import scipy as np
from scipy import linalg, ix_


#===========================================================================
#   MatrixBuilder
#===========================================================================


class MatrixBuilder(object):
    """ Matrix Builder

    Static Members:
        __type_int__ = "Input indices are not int!"

    Instance Members:
        K = global stiffness matrix
        
    Public Methods:
        MatrixBuilder(ndof)
        resize(ndof)
        addValue(idof, jdof, val)
        setValue(idof, jdof, val)
        val = getValue(idof, jdof)
        addBlock(idofs, jdofs, block)
        setBlock(idofs, jdofs, block)
        block = getBlock(idofs, jdofs)
        K = getMatrix()
        print()

    Private Methods:
        __inputIndices(idof, jdof)
        __inputLists(idofs, jdofs)
        
    """
    # Static:
    __type_int__ = "Input indices are not int!"

    # Public:
    def __init__(self,ndof):
        """ Input: ndof = size of square matrix """
        self.K = np.zeros((ndof,ndof),dtype='float64')
    
    def resize(self,ndof):
        """ Input: ndof = new size of square matrix K """
        self.K.resize((ndof,ndof))

    #-----------------------------------------------------------------
    #   add single value
    #-----------------------------------------------------------------

    def setValue(self, idof, jdof, val):
        """
        Input:
            idof = row index
            jdof = col index
            val = value to be set in K[idof,jdof]
        """
        self.K[idof, jdof] = val

    def addValue(self,idof,jdof,val):
        """
        Input:
            idof = row index
            jdof = col index
            val = value to be added to K[idof,jdof]
        """
        self.K[idof,jdof] += val
                
    def getValue(self,idof,jdof):
        """
        Input:
            idof = row index
            jdof = col index
        Output:
            value of K[idof,jdof]  
        """
        return self.K[idof,jdof]

    #-----------------------------------------------------------------
    #   add block
    #-----------------------------------------------------------------
    
    def setBlock(self,idofs,jdofs,block):
        """
        Input:
            idofs = list of row indices
            jdofs = list of col indices
            block = block to be set in K[idofs,jdofs]
        """
        self.K[ix_(idofs,jdofs)] = block
        
    def addBlock(self,idofs,jdofs,block):
        """
        Input:
            idofs = list of row indices
            jdofs = list of col indices
            block = block to be added to K[idofs,jdofs]
        """
        self.K[ix_(idofs,jdofs)] += block
        
    def getBlock(self,idofs,jdofs):
        """
        Input:
            idofs = list of row indices
            jdofs = list of col indices
        Output:
            matrix block K[idofs,jdofs]
        """
        return self.K[ix_(idofs,jdofs)]
    
    def getMatrix(self):
        """ Output: K """
        return self.K
    
    def print(self):
        print("\n",self.K,"\n")
    
    # Private:
    def __inputIndices(self, idof, jdof):
        if isinstance(idof, int) and isinstance(jdof, int):
            pass
        else:
            raise TypeError(self.__type_int__)

    def __inputLists(self, idofs, jdofs):
        # future implementation
        pass


#===========================================================================
#   Other Functions
#===========================================================================


def det(v):
    if isinstance(v,np.ndarray):
        v = linalg.det(v)
    return v

def inv(v):
    if isinstance(v,np.ndarray):
        v = linalg.inv(v)
    else:
        v = 1/v
    return v


#===========================================================================
#   Example
#===========================================================================


if __name__ == "__main__":
    
    mbuild = MatrixBuilder(10)

    # value
    mbuild.setValue(0,0,1)
    mbuild.addValue(0,0,8)
    print("\nValue = ",mbuild.getValue(0,0))

    # block
    block = np.ones((6,6))
    idofs = [3,4,5,6,7,8]
    mbuild.setBlock(idofs, idofs, block)

    idofs = range(1,7)
    mbuild.addBlock(idofs, idofs, block)

    idofs = np.array([2,3,4,5,6,7])
    mbuild.addBlock(idofs, idofs, block)
    print("\nBlock =\n\n",mbuild.getBlock((1,2,3,4),(1,2,3)))

    print("\nSystem Matrix = ")
    mbuild.print()



