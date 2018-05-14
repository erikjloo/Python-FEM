# Import Standard Libraries
import scipy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from copy import deepcopy
from tkinter import Tk, filedialog
from indexed import IndexedOrderedDict

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

        nrow = number of rows (nodes)
        types = list of dof type names
        dofspace = array of dof indices (idofs)
        ndof = last dof index

    Public Methods:
        readMesh(path) - reads gmsh 2.0 file
        readXML(path) - reads .xml file
        plotMesh(rank=2) - plots 2D or 3D
    """

    # Public:

    #-----------------------------------------------------------------------
    #   constructor
    #-----------------------------------------------------------------------

    def __init__(self):

        NodeSet.__init__(self)
        ElementSet.__init__(self)
        self.groups = []
        self.groupNames = IndexedOrderedDict()
        self.ngroups = 0

    #-----------------------------------------------------------------------
    #   initialize
    #-----------------------------------------------------------------------

    def initialize(self, props):
        """ Input:  props = Properties
                    rank = number of dimensions """
        # Read Mesh
        type = props.get("type")
        path = props.get("file")
        rank = props.get("rank")
        doElemGroups = props.get("doElemGroups")
        self.readMesh(type, path, rank, doElemGroups)

        # Initialize DofSpace
        self.rank = rank
        DofSpace.__init__(self, self.nnod, rank)

    #-----------------------------------------------------------------------
    #   readMesh
    #-----------------------------------------------------------------------

    def readMesh(self, type, path, rank, doElemGroups):
        if type == "Gmsh":
            self.readGmsh(path, rank, doElemGroups)
        elif type == "XML":
            self.readXML(path, rank)
        else:
            print("Type can only be Gmsh or XML!")

    #-----------------------------------------------------------------------
    #   readGmsh
    #-----------------------------------------------------------------------

    def readGmsh(self, path=None, rank=3, doElemGroups=False):
        """ Input: path = path_to_file """

        if path is None:
            Tk().withdraw()
            self.path = filedialog.askopenfilename()
        else:
            self.path = path

        fid = open(self.path, "r")

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
                    coord = list(map(float, coord[1:rank+1]))  # coords as float
                    self.addNode(coord)

            #---------------------------------------------------------------
            #   Elements
            #---------------------------------------------------------------

            if line.find('$Elements') == 0:
                data = fid.readline().split()
                nele = int(data[0])

                for iele in range(nele):
                    data = fid.readline().split()
                    ntags = int(data[2])           # number of tags

                    if ntags > 0:
                        Id = int(data[3])       # set group Id
                        if Id not in self.groupNames:
                            groupName = ('Group {}').format(Id)
                            self.groupNames[Id] = groupName
                            print(groupName,"created")
                            self.ngroups += 1
                            self.groups.append([])

                    if doElemGroups is True:
                        Id = self.groupNames.keys().index(Id)
                        self.groups[Id].append(iele)
                    else:
                        self.groups[0].append(iele)

                    connect = list(map(int, data[3+ntags:]))
                    connect = [x-1 for x in connect]
                    self.addElement(connect)

        print(("Mesh read with {} nodes and {} elements").format(
            self.nnod, self.nele))
        fid.close()

    #-----------------------------------------------------------------------
    #   readXML
    #-----------------------------------------------------------------------

    def readXML(self, path=None, rank=3):
        """ Input: path = path_to_file """
        if path is None:
            Tk().withdraw()
            self.path = filedialog.askopenfilename()
        else:
            self.path = path

        with open(path, 'r') as file:

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

                        coord = list(map(float, data[1:rank+1]))
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

    #-----------------------------------------------------------------------
    #   plotDeformed
    #-----------------------------------------------------------------------

    def plotDeformed(self, disp, scale, rank=2):
        """ Input:  disp = displacement vector
                    scale = 
                    rank = number of dimensions """
        
        # Craft deformed coordinates
        deformed = self.getCoords()
        for inod in range(len(self.coords)):
            idofs = self.getDofIndices(inod)
            x = self.getCoords(inod)
            u = np.array(disp[idofs])*scale
            deformed[inod,:] = u+x
        
        # Create figure
        if rank == 1 or rank == 2:
            fig = plt.figure(figsize=(6, 6))
            ax = fig.add_subplot(111)
        elif rank == 3:
            fig = plt.figure(figsize=(6, 6))
            ax = fig.add_subplot(111, projection='3d')

        # Plot deformed mesh
        for connect in self.connectivity:
            coords = deformed[np.ix_(connect), :][0]
            ax.plot(coords[:, 0], coords[:, 1], linewidth=0.5, color='k')

        plt.show()

    #-----------------------------------------------------------------------
    #   updateGeometry
    #-----------------------------------------------------------------------

    def updateGeometry(self, disp):
        """ Input:  disp = displacement vector """
        for inod in range(len(self.coords)):
            idofs = self.getDofIndices(inod)
            x = [a+b for a,b in zip(self.coords[inod],disp[idofs])]
            self.coords[inod] = x


#===========================================================================
#   Example
#===========================================================================


if __name__ == '__main__':

    from properties import Properties
    
    mesh = Mesh()
    mesh.readXML("Examples/square.xml")
    mesh.plotMesh(rank=2)

    file = "Examples/semicircle.pro"
    props = Properties()
    props.parseFile(file)

    mesh = Mesh()
    mesh.initialize(props.getProps("input.mesh"))
    mesh.plotMesh(rank=2)