# Import Standard Libraries
import scipy as np
import matplotlib.pyplot as plt
from pprint import PrettyPrinter
from scipy.sparse import dok_matrix
from scipy.linalg import inv, det, norm, svd, eigvals, cholesky, LinAlgError

#===========================================================================
#   MatrixBuilder
#===========================================================================


class MatrixBuilder(object):
    """ Dok Sparse Matrix Builder

    Static Members:
        __type_int__ = "Input idof or jdof is not int!"

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
        K = getDenseMatrix()
        K = getMatrix()
        print()
        plot()

    """

    # Static:
    __type_int__ = "Input idof or jdof is not int!"

    # Public:

    #-----------------------------------------------------------------
    #   constructor
    #-----------------------------------------------------------------
    
    def __init__(self,ndof):
        """ Input: ndof = size of square matrix """
        self.hbw = 0
        self.K = dok_matrix((ndof,ndof),dtype="float64")
    
    def resize(self,ndof):
        """ Input: ndof = new size of square matrix K """
        self.K.resize((ndof,ndof))

    #-----------------------------------------------------------------
    #   add single value
    #-----------------------------------------------------------------

    def setValue(self, idof, jdof, val):
        """ Input:  idof = row index
                    jdof = col index
                    val = value to be set in K[idof,jdof] """
        if isinstance(idof, int) and isinstance(jdof, int):
            self.K[idof, jdof] = val
        else:
            raise TypeError(self.__type_int__)

    def addValue(self,idof,jdof,val):
        """ Input:  idof = row index
                    jdof = col index
                    val = value to be added to K[idof,jdof] """
        if isinstance(idof, int) and isinstance(jdof, int):
            self.K[idof,jdof] += val
        else:
            raise TypeError(self.__type_int__)

    def getValue(self,idof,jdof):
        """ Input:  idof = row index
                    jdof = col index
            Output: value of K[idof,jdof] """
        return self.K[idof,jdof]

    #-----------------------------------------------------------------
    #   add block
    #-----------------------------------------------------------------
    
    def setBlock(self,idofs,jdofs,block):
        """ Input:  idofs = list of row indices
                    jdofs = list of col indices
                    block = block to be set in K[idofs,jdofs] """
        self.K[np.ix_(idofs,jdofs)] = block
        
    def addBlock(self,idofs,jdofs,block):
        """ Input:  idofs = list of row indices
                    jdofs = list of col indices
                    block = block to be added to K[idofs,jdofs] """
        self.K[np.ix_(idofs,jdofs)] += block
        
    def getBlock(self,idofs,jdofs):
        """ Input:  idofs = list of row indices
                    jdofs = list of col indices
            Output: matrix block K[idofs,jdofs] """
        return self.K[np.ix_(idofs,jdofs)].todense()

    #-----------------------------------------------------------------
    #   get matrix
    #-----------------------------------------------------------------

    def getDenseMatrix(self):
        """ Output: K """
        return self.K.todense()

    def getMatrix(self):
        """ Output: K """
        return self.K

    #-----------------------------------------------------------------
    #   print and plot
    #-----------------------------------------------------------------

    def print(self):
        pp = PrettyPrinter(indent=1, width=120, compact=True)
        pp.pprint(self.K.todense()) 
    
    def plot(self):
        plt.spy(self.K)
        plt.show()

#===========================================================================
#   Other Functions
#===========================================================================


def determinant(v):
    if isinstance(v,np.ndarray):
        return det(v)
    else:
        return np.array(v)

def inverse(v):
    if isinstance(v,np.ndarray):
        return inv(v)
    else:
        return np.array(1/v)

def gram_schmidt(u1, u2, u3=None):
    if u3 is None:
        v1 = u1/norm(u1)
        v2 = u2 - np.dot(u2,v1)/np.dot(v1,v1)*v1
        v2 = v2/norm(v2)
        return v1, v2
    else:
        v1 = u1/norm(u1)
        v2 = u2 - np.dot(u2,v1)/np.dot(v1,v1)*v1
        v3 = u3 - np.dot(u3,v1)/np.dot(v1,v1)*v1 - \
                - np.dot(u3, v2)/np.dot(v2, v2)*v2
        v2 = v2/norm(v2)
        v3 = v3/norm(v3)
        return v1, v2, v3


#===========================================================================
#   Example
#===========================================================================


if __name__ == "__main__":
    
    mbuild = MatrixBuilder(10)

    # value
    mbuild.setValue(0,0,1)
    mbuild.addValue(0,0,8)
    print("\nValue = ",mbuild.getValue(0,0))

    # Add block 1
    block = np.ones((6,6))
    idofs = [3,4,5,6,7,8]
    mbuild.setBlock(idofs, idofs, block)

    # Add block 2
    idofs = range(1,7)
    mbuild.addBlock(idofs, idofs, block)

    # Add block 3
    idofs = np.array([2,3,4,5,6,7])
    mbuild.addBlock(idofs, idofs, block)
    print("\nBlock =\n\n",mbuild.getBlock((1,2,3,4),(1,2,3)))

    print("\nSystem Matrix = ")
    mbuild.print()

    print("\nSystem Matrix = ")
    mbuild.resize(5)
    mbuild.print()

    mbuild.plot()



