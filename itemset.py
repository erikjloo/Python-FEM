# Import Standard Libraries
import scipy as np
from itertools import chain


#===========================================================================
#   ItemSet
#===========================================================================


class ItemSet(object):
    """ Item Set

    Static Members:
        __type__ = "Input is not list or array!"
        __type_int__ = "Input is not int!"
        __type_int_list__ = "Input is not int or list!"
    """

    # Static:
    __type__ = "Input is not list or array!"
    __type_int__ = "Input is not int!"
    __type_int_list__ = "Input is not int or list!"


#===========================================================================
#   NodeSet
#===========================================================================


class NodeSet(ItemSet):
    """ Node Set

    Static Members:
        __type__ = "Input is not list or array!"
        __type_int__ = "Input is not int!"
        __type_int_list__ = "Input is not int or list!"

    Instance Members:
        coords = list of nodal coordinates
        nnod = number of nodes

    Public Methods:
        NodeSet()
        inod = addNode(coord)
        inodes = addNodes(coords)
        setNode(inod, coord)
        eraseNode(inod)
        eraseNodes(inodes)
        nnod = nodeCount()
        coord[s] = getCoords(inod[es])
    """

    # Public:
    def __init__(self):
        self.coords = []
        self.nnod = 0

    def addNode(self, coord):
        """ Input: coord = coordinates of new node 
            Output: inod = node index """
        if isinstance(coord, list):
            self.coords.append(coord)
        elif isinstance(coord, np.ndarray):
            self.coords.append(coord.tolist())
        else:
            raise TypeError(self.__type__)
        self.nnod += 1
        return self.nnod - 1

    def addNodes(self, coords):
        """ Input: coords = list of coordinates of new nodes 
            Output: inodes = list of node indices """
        inodes = []
        if isinstance(coords, list):
            for coord in coords:
                inodes.append(self.addNode(coord))
        elif isinstance(coords, np.ndarray):
            for coord in coords.tolist():
                inodes.append(self.addNode(coord))
        else:
            raise TypeError(self.__type__)
        return inodes

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
        if isinstance(inod, int):
            del self.coords[inod]
            self.nnod -= 1
        else:
            raise TypeError(self.__type_int__)

    def eraseNodes(self, inodes):
        """ Input: inodes = (list of) indices of nodes to be erased """
        if isinstance(inodes, (list, tuple, range, np.ndarray)):
            for inod in sorted(inodes, reverse=True):
                self.eraseNode(inod)
        elif isinstance(inodes, int):
            self.eraseNode(inodes)
        else:
            raise TypeError(self.__type_int_list__)

    def nodeCount(self):
        """ Output: number of nodes """
        return self.nnod

    def getCoords(self, inodes=None):
        """ Input: inodes = (list of) node indices
            Output: coordinates of inodes (if given) """
        if inodes is None:
            return np.array(self.coords)
        elif isinstance(inodes, (list, tuple, range, np.ndarray)):
            return np.array([self.coords[inod] for inod in inodes])
        elif isinstance(inodes, int):
            return np.array(self.coords[inodes])
        else:
            raise TypeError(self.__type_int_list__)


#===========================================================================
#   ElementSet
#===========================================================================


class ElementSet(ItemSet):
    """ Element set

    Static Members:
        __type__ = "Input is not list or array!"
        __type_int__ = "Input is not int!"
        __type_int_list__ = "Input is not int or list!"

    Instance Members:
        connectivity = list of element connectivities
        nele = number of elements

    Public Methods:
        ElementSet()
        iele = addElement(connect)
        ielements = addElements(connect)
        setElement(iele, connect)
        eraseElement(iele)
        eraseElements(ielements)
        nele = elemCount()
        connect[ivity] = getNodes(iele[ments])
        inodes = getNodeIndices(iele[ments])
    """

    # Public:
    def __init__(self):
        self.connectivity = []
        self.nele = 0

    def addElement(self, connect):
        """ Input: connect = node indices of new element 
            Output: iele = element index """
        if isinstance(connect, list):
            self.connectivity.append(connect)
        elif isinstance(connect, np.ndarray):
            self.connectivity.append(connect.tolist())
        else:
            raise TypeError(self.__type__)
        self.nele += 1
        return self.nele - 1

    def addElements(self, connectivity):
        """ Input: connectivity = list of node indices of new elements 
            Output: ielements = list of element indices """
        ielements = []
        if isinstance(connectivity, list):
            for connect in connectivity:
                ielements.append(self.addElement(connect))
        elif isinstance(connectivity, np.ndarray):
            for connect in connectivity.tolist():
                ielements.append(self.addElement(connect))
        else:
            raise TypeError(self.__type__)
        return ielements

    def setElement(self, iele, connect):
        """ Input: iele = element index, connect = element node indices """
        if isinstance(connect, list):
            self.connectivity[iele] = connect
        elif isinstance(connect, np.ndarray):
            self.connectivity[iele] = connect.tolist()
        else:
            raise TypeError(self.__type__)

    def eraseElement(self, iele):
        """ Input: iele = index of element to be erased """
        if isinstance(inod, int):
            del self.connectivity[iele]
            self.nele -= 1
        else:
            raise TypeError(self.__type_int__)

    def eraseElements(self, ielements):
        """ Input: ieles = (list of) indices of elements to be erased """
        if isinstance(ielements, (list, tuple, range, np.ndarray)):
            for iele in sorted(ielements, reverse=True):
                self.eraseElement(iele)
        elif isinstance(ielements, int):
            self.eraseElement(ielements)
        else:
            raise TypeError(self.__type_int_list__)

    def elemCount(self):
        """ Output: number of elements """
        return self.nele

    def getNodes(self, ielements=None):
        """ Input: ielements = ielements = (list of) element indices
            Output: element connectivity vector(s) of ielements (if given) """
        if ielements is None:
            return self.connectivity
        elif isinstance(ielements, (list, tuple, range, np.ndarray)):
            connectivity = []
            for iele in ielements:
                connectivity.append(self.connectivity[iele])
            return connectivity
        elif isinstance(ielements, int):
            return self.connectivity[ielements]
        else:
            raise TypeError(self.__type_int_list__)
    
    def getNodeIndices(self, ielements=None):
        """ Input: ielements = (list of) element indices
            Output: inodes = node indices of given element indices """
        connectivity = self.getNodes(ielements)
        return list(set(chain.from_iterable(connectivity)))


#===========================================================================
#   Example
#===========================================================================


if __name__ == '__main__':

    nodes = NodeSet()
    inod = nodes.addNode([0, 0, 0])  # works with list and array
    print("\nNode ", inod, " created.\n")
    inodes = nodes.addNodes(
        np.array([[1, 0, 0], [1, 1, 0], [0.5, 2, 0], [0, 1, 0]]))
    print("Nodes ", inodes, " created.\n")
    nodes.setNode(3, [0.5, 3, 0, ])
    print("Nodes :\n", nodes.getCoords())
    print("Node count :", nodes.nodeCount(), "\n")

    nodes.eraseNode(4)  # eraseNodes(4) works too!
    nodes.eraseNodes(range(2, 4))  # works with list, tuple, range & array
    inodes = nodes.addNodes([[1, 1, 0], [0.5, 2, 0], [0, 1, 0]])
    print("Nodes ", inodes, " created again.\n")
    print("Nodes :\n", nodes.getCoords(np.array([0, 1, 2, 3, 4])))
    print("Node count :", nodes.nodeCount(), "\n")

    elems = ElementSet()
    iele = elems.addElement([0, 1])  # works with list and array
    print("\nElement ", iele, " created.\n")
    ielements = elems.addElements(np.array([[1, 2], [2, 4], [3, 4], [4, 0]]))
    print("Elements ", ielements, " created.\n")
    elems.setElement(2, [2, 3])
    print("Elements :\n", elems.getNodes())
    print("Elem count :", elems.elemCount(), "\n")

    elems.eraseElement(4) # works with list, tuple, range & array
    elems.eraseElements(np.array([2, 3]))
    ielements = elems.addElements([[2, 3], [3, 4], [4, 0]])
    print("\nElements ", ielements, " created again.\n")
    print("Elements :\n", elems.getNodes(range(5)))
    print("Elem count :", elems.elemCount(), "\n")

    inodes = elems.getNodeIndices(ielements)
    print("\n Node indices :\n", inodes)
