# Import Standard Libraries
import scipy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from tkinter import Tk, filedialog

#Import Local Libraries
from itemset import NodeSet, ElementSet
from dofspace import DofSpace


#===========================================================================
#   Mesh
#===========================================================================


class Mesh(NodeSet, ElementSet, DofSpace):
    """ Mesh class

    Static Members:
        __type__ = "Input is not list or array!"

    Instance Members:
        coords = list of nodal coordinates
        nnod = last node index + 1

        connectivity = list of element connectivities
        nele = last element index + 1

        groups = list of elements in each physical group
        groupNames = names of physical groups
        ngroups = number of physical groups

        dofspace = array of dof indices (idofs)
        types = list of dof type names
        idof = last dof index = ndof - 1
        ndof = number of dofs

    Public Methods:
        readMesh(__path__) - reads gmsh 2.0 file
        readXML(__path__) - reads .xml file
        plotMesh(rank=2) - plots 2D or 3D

    """
    # Public:

    #=======================================================================
    #   constructor
    #=======================================================================

    def __init__(self):

        NodeSet.__init__(self)
        ElementSet.__init__(self)
        self.groups = []
        self.groupNames = {}
        self.ngroups = 0

    #-----------------------------------------------------------------------
    #   initialize
    #-----------------------------------------------------------------------

    def initialize(self, props, rank=1):

        # Read Mesh
        type = props.get("Input.mesh.type")
        path = props.get("Input.mesh.file")
        self.readMesh(type, path)

        # Initialize DofSpace
        self.rank = rank
        DofSpace.__init__(self, self.nnod, rank)

    #-----------------------------------------------------------------------
    #   readMesh
    #-----------------------------------------------------------------------

    def readMesh(self, type, __path__):
        if type == "Gmsh":
            self.readGmsh(__path__)
        elif type == "XML":
            self.readXML(__path__)
        else:
            print("Type can only be Gmsh or XML!")

    #-----------------------------------------------------------------------
    #   readGmsh
    #-----------------------------------------------------------------------

    def readGmsh(self, __path__=None):
        """ Input: __path__ = path_to_file """

        if __path__ is None:
            Tk().withdraw()
            self.__path__ = filedialog.askopenfilename()
        else:
            self.__path__ = __path__

        fid = open(self.__path__, "r")

        line = "start"

        while line:
            line = fid.readline()

            #---------------------------------------------------------------
            #   Groups
            #---------------------------------------------------------------

            if line.find('$PhysicalNames') == 0:
                data = fid.readline().split()
                self.ngroups = int(data[0])
                self.groups = [[] for _ in range(self.ngroups)]

                for _ in range(self.ngroups):
                    line = fid.readline()
                    newkey = int(line.split()[0])
                    qstart = line.find('"')+1
                    qend = line.find('"', -1, 0)-1
                    self.groupNames[newkey] = line[qstart:qend]

            #---------------------------------------------------------------
            #   Nodes
            #---------------------------------------------------------------

            if line.find('$Nodes') == 0:
                data = fid.readline().split()
                nnod = int(data[0])

                for _ in range(nnod):
                    coord = fid.readline().split()      # coords as str
                    coord = list(map(float, coord[1:]))  # coords as float
                    self.addNode(coord)

            #---------------------------------------------------------------
            #   Elements
            #---------------------------------------------------------------

            if line.find('$Elements') == 0:
                data = fid.readline().split()
                nele = int(data[0])

                for iele in range(nele):
                    data = fid.readline().split()
                    etype = int(data[1])           # element type
                    ntags = int(data[2])           # number of tags

                    if ntags > 0:
                        Id = int(data[3])       # set group Id
                        if Id not in self.groupNames:
                            self.groupNames[Id] = ('Group {}').format(Id)
                            self.ngroups += 1
                            self.groups.append([])

                    self.groups[Id].append(iele)

                    connect = list(map(int, data[3+ntags:]))
                    connect = [x-1 for x in connect]
                    self.addElement(connect)

        print(("Mesh read with {} nodes and {} elements").format(
            self.nnod, self.nele))
        fid.close()

    #-----------------------------------------------------------------------
    #   readXML
    #-----------------------------------------------------------------------

    def readXML(self, __path__=None):
        """ Input: __path__ = path_to_file """
        if __path__ is None:
            Tk().withdraw()
            self.__path__ = filedialog.askopenfilename()
        else:
            self.__path__ = __path__

        with open(__path__, 'r') as file:

            flag_n = False
            flag_e = False

            for line in file:

                if line.startswith("<Nodes>"):
                    flag_n = True
                elif line.startswith("</Nodes>"):
                    flag_n = False
                elif line.startswith("<Elements>"):
                    flag_e = True
                elif line.startswith("</Elements>"):
                    flag_e = False

                data = line.split(';')
                data = data[0].split()

                if len(data) > 0 and data[0].isdigit():

                    #-------------------------------------------------------
                    #   Nodes
                    #-------------------------------------------------------

                    if flag_n is True:

                        coord = list(map(float, data[1:]))
                        self.addNode(coord)

                    #-------------------------------------------------------
                    #   Elements
                    #-------------------------------------------------------

                    if flag_e is True:

                        connect = list(map(int, data[1:]))
                        connect = [x-1 for x in connect]
                        self.addElement(connect)

        print(("Mesh read with {} nodes and {} elements").format(
            self.nnod, self.nele))

    #-----------------------------------------------------------------------
    #   plotMesh
    #-----------------------------------------------------------------------

    def plotMesh(self, rank=2):
        """ Input: rank = number of dimensions """
        # Determine whether 2D or 3D plot
        if rank == 1 or rank == 2:
            fig = plt.figure(figsize=(6, 6))
            ax = fig.add_subplot(111)
        elif rank == 3:
            fig = plt.figure(figsize=(6, 6))
            ax = fig.add_subplot(111, projection='3d')

        # Plot Mesh
        for connect in self.connectivity:
            coords = self.getCoords(connect)
            ax.plot(coords[:, 0], coords[:, 1], linewidth=0.5, color='k')

        plt.show()


#===========================================================================
#   Example
#===========================================================================


if __name__ == '__main__':

    from properties import Properties
    
    mesh = Mesh()
    mesh.readXML("Examples/square.xml")
    mesh.plotMesh(rank=2)

    file = "Examples/square.pro"
    props = Properties()
    props.parseFile(file)

    mesh = Mesh()
    mesh.initialize(props)
    mesh.printDofSpace()
    mesh.plotMesh()
