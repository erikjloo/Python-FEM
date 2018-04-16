# Import Standard Libraries
import scipy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from tkinter import Tk, filedialog

#Import Local Libraries
from itemset import NodeSet, ElementSet


#===========================================================================
#   Mesh
#===========================================================================


class Mesh(NodeSet, ElementSet):
    """ Mesh class

    Static Members:
        __type__ = "Input is not list or array!"

    Instance Members:
        coords = list of nodal coordinates
        inod = last node index = nnod - 1
        nnod = number of nodes

        connectivity = list of element connectivities
        iele = last element index = nele - 1
        nele = number of elements

        props = list of element type and physical group of each element
        Phys = names of physical groups
        nPhys = number of physical groups

    Public Methods:
        readMesh(__path__)
        plotMesh()

    """
    # Public:

    def __init__(self):

        self.coords = []
        self.inod = -1
        self.nnod = 0

        self.connectivity = []
        self.iele = -1
        self.nele = 0

        self.props = []
        self.Phys = {}
        self.nPhys = 0

    #=============================================================
    #   readMesh
    #=============================================================

    def readMesh(self, __path__=None):
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

            #-----------------------------------------------------
            #   Physical Names
            #-----------------------------------------------------

            if line.find('$PhysicalNames') == 0:
                data = fid.readline().split()
                self.nPhys = int(data[0])

                for _ in range(self.nPhys):
                    line = fid.readline()
                    newkey = int(line.split()[0])
                    qstart = line.find('"')+1
                    qend = line.find('"', -1, 0)-1
                    self.Phys[newkey] = line[qstart:qend]

                if fid.readline().find('$EndPhysicalNames') != 0:
                    raise ValueError('expecting EndPhysicalNames')

            #-----------------------------------------------------
            #   Nodes
            #-----------------------------------------------------

            if line.find('$Nodes') == 0:
                nnod = fid.readline().split()
                self.nnod = int(nnod[0])
                self.inod += self.nnod

                for _ in range(self.nnod):
                    coord = fid.readline().split()      # coords as str
                    coord = list(map(float, coord[1:]))  # coords as float
                    self.coords.append(coord)

                if fid.readline().find('$EndNodes') != 0:
                    raise ValueError('expecting EndNodes')

            #-----------------------------------------------------
            #   Elements
            #-----------------------------------------------------

            if line.find('$Elements') == 0:
                nele = fid.readline().split()
                self.nele = int(nele[0])
                self.iele += self.nele

                for _ in range(self.nele):
                    data = fid.readline().split()
                    etype = int(data[1])           # element type
                    ntags = int(data[2])           # number of tags

                    if ntags > 0:
                        physid = int(data[3])       # set physical id
                        if physid not in self.Phys:
                            self.Phys[physid] = (
                                'Physical Entity {}').format(physid)
                            self.nPhys += 1

                    self.props.append([etype, physid])

                    connect = list(map(int, data[3+ntags:]))
                    connect = [x-1 for x in connect]
                    self.connectivity.append(connect)

                line = fid.readline()
                if line.find('$EndElements') != 0:
                    raise ValueError('expecting EndElements')
        fid.close()

    #=============================================================
    #   readXML
    #=============================================================

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
                    self.nnod = self.inod + 1
                elif line.startswith("<Elements>"):
                    flag_e = True
                elif line.startswith("</Elements>"):
                    flag_e = False
                    self.nele = self.iele + 1

                data = line.split(';')
                data = data[0].split()

                if len(data) > 0 and data[0].isdigit():

                    #-----------------------------------------------------
                    #   Nodes
                    #-----------------------------------------------------

                    if flag_n is True:

                        coord = list(map(float, data[1:]))
                        self.coords.append(coord)
                        self.inod += 1

                    #-----------------------------------------------------
                    #   Elements
                    #-----------------------------------------------------

                    if flag_e is True:

                        connect = list(map(int, data[1:]))
                        connect = [x-1 for x in connect]
                        self.connectivity.append(connect)
                        self.iele += 1

    #=============================================================
    #   plotMesh
    #=============================================================

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

    mesh = Mesh()
    mesh.readXML("example.xml")

    coords = mesh.getCoords()
    print(coords)

    connect = mesh.getNodes()
    print(connect)

    mesh.plotMesh(rank=2)

    mesh1 = Mesh()
    mesh1.readMesh("square.msh")

    coords = mesh1.getCoords(range(10))
    print(coords)

    connect = mesh1.getNodes(range(10))
    print(connect)

    mesh1.plotMesh()
