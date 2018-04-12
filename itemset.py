# Import Standard Libraries
import scipy as np
import warnings
from scipy import ix_


#===========================================================================
#   NodeSet
#===========================================================================


class NodeSet(object):
    """ Node Set

    Static Members:
        __type__ = "Input is not list or array!"
    Instance Members:
        coords = list of nodal coordinates
    Public Methods:
        NodeSet()
        addNode(coord)
        addNodes(coords)
        setNode(coord)
        eraseNode(inod)
        eraseNodes(inodes)
        nnod = nodeCount()
        coord[s] = getCoords(inod[es])
    """
    # Static:
    __type__ = "Input is not list or array!"
    
    # Public:
    def __init__(self):
        self.coords = []
        
    def addNode(self, coord):
        """ Input: coord = coordinates of new node """
        if isinstance(coord, list):
            self.coords.append(coord)
        elif isinstance(coord, np.ndarray):
            self.coords.append(coord.tolist())
        else:
            raise TypeError(self.__type__)

    def addNodes(self, coords):
        """ Input: coords = list of coordinates of new nodes """
        if isinstance(coords, list):
            for coord in coords:
                self.addNode(coord)
        elif isinstance(coords, np.ndarray):
            for coord in coords.tolist():
                self.addNode(coord)
        else:
            raise TypeError(self.__type__)
                
    def setNode(self, inod, coord):
        """ Input: inod = node index, coord = node coodinates """
        if isinstance(coord, list):
            self.coords[inod] = coord
        elif isinstance(coord, np.ndarray):
            self.coords[inod] = coord.tolist()
        else:
            raise TypeError(self.__type__)
        
    def eraseNode(self, inod):
        """ Input: inod = index of node to be erased """
        del self.coords[inod]
        
    def eraseNodes(self, inodes):
        """ Input: inodes = indices of nodes to be erased """
        if isinstance(inodes, (list,tuple,range,np.ndarray)):
            for inod in sorted(inodes, reverse=True):
                self.eraseNode(inod)
        else:
            self.eraseNode(inodes)

    def nodeCount(self):
        """ Output: number of nodes """
        return len(self.coords)

    def getCoords(self, inodes=None):
        """
        Input:
            inodes = node index or node indices
        Output:
            coordinates of inodes (if given)
        """
        coords = np.array(self.coords)
        if inodes is None:
            return coords
        elif isinstance(inodes, (list,tuple,range,np.ndarray)):
            return coords[ix_(inodes),:][0]
        else:
            return coords[inodes]


#===========================================================================
#   ElementSet
#===========================================================================


class ElementSet(object):
    """ Element set

    Static Members:
        __type__ = "Input is not list or array!"
    Instance Members:
        connectivity = list of element connectivities
    Public Methods:
        ElementSet()
        addElement(connect)
        addElements(connect)
        setElement(connect)
        eraseElement(iele)
        eraseElements(ielements)
        nele = elemCount()
        connect[ivity] = getNodes(iele[ments])
    """
    # Static:
    __type__ = "Input is not list or array!"
    
    # Public:
    def __init__(self):
        self.connectivity = []
        
    def addElement(self, connect):
        """ Input: connect = node indices of new element """
        if isinstance(connect, list):
            self.connectivity.append(connect)
        elif isinstance(connect, np.ndarray):
            self.connectivity.append(connect.tolist())
        else:
            raise TypeError(self.__type__)

    def addElements(self, connectivity):
        """ Input: connectivity = list of node indices of new elements """
        if isinstance(connectivity, list):
            for connect in connectivity:
                self.addElement(connect)
        elif isinstance(connectivity, np.ndarray):
            for connect in connectivity.tolist():
                self.addElement(connect)
        else:
            raise TypeError(self.__type__)
            
    def setElement(self, iele, connect):
        """ Input: iele = element index, connect = element nodes indices """
        if isinstance(connect, list):
            self.connectivity[iele] = connect
        elif isinstance(connect, np.ndarray):
            self.connectivity[iele] = connect.tolist()
        else:
            raise TypeError(self.__type__)
        
    def eraseElement(self, iele):
        """ Input: iele = index of element to be erased """
        del self.connectivity[iele]
    
    def eraseElements(self, ielements):
        """ Input: ieles = indices of elements to be erased """
        if isinstance(ielements, (list,tuple,range,np.ndarray)):
            for iele in sorted(ielements, reverse=True):
                self.eraseElement(iele)
        else:
            self.eraseElement(ielements)
            
    def elemCount(self):
        """ Output: number of elements """
        return len(self.connectivity)
    
    def getNodes(self, ielements=None):
        """
        Input:
            ielements = element index or element indices
        Output:
            element connectivity vector(s) of ielements (if given)
        """
        # connectivity = np.array(self.connectivity)
        if ielements is None:
            return self.connectivity
        elif isinstance(ielements, (list,tuple,range,np.ndarray)):
            connectivity = []
            for iele in ielements:
                connectivity.append(self.connectivity[iele])
            return connectivity
        else:
            return self.connectivity[ielements]


#===========================================================================
#   Example
#===========================================================================


if __name__ == '__main__':

    nodes = NodeSet()
    nodes.addNode([0,0,0]) # works with list and array
    nodes.addNodes(np.array([[1,0,0],[1,1,0],[0.5,2,0],[0,1,0]]))
    nodes.setNode(3,[0.5,3,0,])
    print("Nodes :\n",nodes.getCoords())
    print("Node count :", nodes.nodeCount())
    
    nodes.eraseNode(4) # eraseNodes(4) works too!
    nodes.eraseNodes(range(2,4)) # works with list, tuple, range & array
    nodes.addNodes([[1,1,0],[0.5,2,0],[0,1,0]])
    print("\nNodes :\n",nodes.getCoords(np.array([0,1,2,3,4])))
    print("Node count :", nodes.nodeCount())

    elems = ElementSet()
    elems.addElement([0,1]) # works with list and array
    elems.addElements(np.array([[1,2],[2,4],[3,4],[4,0]]))
    elems.setElement(2,[2,3])
    print("\nElements :\n",elems.getNodes())
    print("Elem count :",elems.elemCount())

    elems.eraseElement(4)
    elems.eraseElements(np.array([2,3])) # works with list, tuple, range & array
    elems.addElements([[2,3],[3,4],[4,0]])
    print("\nElements :\n",elems.getNodes(range(5)))
    print("Elem count :",elems.elemCount())
