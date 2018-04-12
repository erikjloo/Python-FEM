# Import Standard Libraries
import scipy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
from mpl_toolkits.mplot3d import proj3d
from tkinter import Tk, filedialog

#Import Local Libraries
from itemset import NodeSet, ElementSet


def orthogonal_plot():
    fig = plt.figure(figsize=(6, 6))
    return fig.add_subplot(111, projection='3d')

    
#===========================================================================
#   Mesh
#===========================================================================


class Mesh(NodeSet,ElementSet):
    """ Mesh class

    Instance Members:
        coords
        connectivity
        props
        Phys
        nnodes
        nele
        nprops
    Public Methods:
        readMesh(__path__)
        plotMesh()
    """
    
    # Public:
    def __init__(self):

        self.coords = []
        self.connectivity = []
        self.props = []
        self.Phys = {}

        self.nnodes = 0
        self.nele = 0
        self.nprops = 0

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
                self.nprops = int(data[0])

                for _ in range(self.nprops):
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
                nnodes = fid.readline().split()
                self.nnodes = int(nnodes[0])

                for _ in range(self.nnodes):
                    coord = fid.readline().split()      # coords as str
                    coord = list(map(float, coord[1:])) # coords as float
                    self.coords.append(coord)

                if fid.readline().find('$EndNodes') != 0:
                    raise ValueError('expecting EndNodes')
                
            #-----------------------------------------------------
            #   Elements
            #-----------------------------------------------------
            
            if line.find('$Elements') == 0:
                nele = fid.readline().split()
                self.nele = int(nele[0])

                for _ in range(self.nele):
                    data = fid.readline().split()
                    etype = int(data[1])           # element type
                    ntags = int(data[2])           # number of tags
                    
                    if ntags > 0:                   
                        physid = int(data[3])       # set physical id
                        if physid not in self.Phys:
                            self.Phys[physid] = ('Physical Entity {}').format(physid)
                            self.nprops += 1

                    self.props.append([etype,physid])
                    
                    connect = list(map(int, data[3+ntags:]))
                    connect = [x-1 for x in connect]
                    self.connectivity.append(connect)

                line = fid.readline()
                if line.find('$EndElements') != 0:
                    raise ValueError('expecting EndElements')
        fid.close()

    def plotMesh(self):
        
        #ax = orthogonal_plot()
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111)
        
        for connect in self.connectivity:
            coords = self.getCoords(connect)
            ax.plot(coords[:,0],coords[:,1], linewidth = 0.5, color='k')
        
        plt.show()
    
#===========================================================================
#   Example
#===========================================================================


if __name__ == '__main__':

    mesh = Mesh()
    mesh.readMesh("rve.msh")
    mesh.plotMesh()
    
    coords = mesh.getCoords(range(10))
    print(coords)

    connect = mesh.getNodes(range(10))
    print(connect)
