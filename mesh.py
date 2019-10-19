# Import Standard Libraries
import re
import logging
import scipy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from tkinter import Tk, filedialog
from indexed import IndexedOrderedDict
# from collections import OrderedDict

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
        __type_int__ = "Input is not int!"
        __type_str__ = "Input is not str!"
        __type_int_list__ = "Input is not int or list!"
        __type_str_list__ = "Input is not str or list!"
        __type_dof__ = "Input inod is not int or dof is not str!"
        __renumber__ = "Erasing dofs: Dof numbers will be renumbered!"
        
    Instance Members:
        coords = list of nodal coordinates
        nnod = number of nodes

        connectivity = list of element connectivities
        nele = number of elements

        nrow = number of rows (nodes)
        types = list of dof names
        dofspace = array of dof indices (idofs)
        ndof = number of degrees of freedom

        groups = list of elements in each physical group
        groupNames = names of physical groups
        ngroups = number of physical groups

        path = file path
        type = file type
        rank = number of dimensions
        doElemGroups = bool
        
    Public Methods:
        Mesh(conf, props)
        readMesh(self, type, path, rank, doElemGroups)
        readGmsh(self, path, rank, doElemGroups)
        readXML(self, path, rank)
        plotMesh(rank)
        plotDeformed(self, disp, scale, rank)
        updateGeometry(self, disp)
    """

    # Public:

    #-----------------------------------------------------------------------
    #   constructor
    #-----------------------------------------------------------------------

    def __init__(self, conf=None, props=None):
        """ Input:  conf = output properties
                    props = input properties """

        NodeSet.__init__(self)
        ElementSet.__init__(self)
        self.groups = []
        self.groupNames = IndexedOrderedDict()
        self.ngroups = 0

        if conf and props:
            # Get Props
            myProps = props.getProps("mesh")
            myConf = conf.makeProps("mesh")

            self.path = myProps.get("file")
            self.type = myProps.get("type", "Gmsh")
            self.rank = myProps.get("rank", 2)
            self.doElemGroups = myProps.get("doElemGroups", False)

            myConf.set("file", self.path)
            myConf.set("type", self.type)
            myConf.set("rank", self.rank)
            myConf.set("doElemGroups", self.doElemGroups)

            # Read Mesh
            self.readMesh(self.type, self.path, self.rank, self.doElemGroups)

            # Initialize DofSpace
            DofSpace.__init__(self, self.nnod, self.rank)

    #-----------------------------------------------------------------------
    #   readMesh
    #-----------------------------------------------------------------------

    def readMesh(self, type, path, rank, doElemGroups):
        if type == "Gmsh":
            self.readGmsh(path, rank, doElemGroups)
        elif type == "XML":
            self.readXML(path, rank, doElemGroups)
        else:
            raise ValueError("type can only be Gmsh or XML!")

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
                    key = int(line.split()[1])
                    qstart = line.find('"')+1
                    qend = line.find('"', -1, 0)-1
                    self.groupNames[key] = line[qstart:qend]

            #---------------------------------------------------------------
            #   Nodes
            #---------------------------------------------------------------

            if line.find('$Nodes') == 0:
                data = fid.readline().split()
                nnod = int(data[0])

                for _ in range(nnod):
                    data = fid.readline().split()
                    coord = [float(x) for x in data[1:rank+1]]
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
                        key = int(data[3])
                        if key not in self.groupNames:
                            groupName = 'Group {}'.format(key)
                            self.groupNames[key] = groupName
                            self.ngroups += 1
                            self.groups.append([])

                    if doElemGroups is True:
                        idx = self.groupNames.keys().index(key)
                        self.groups[idx].append(iele)
                    else:
                        self.groups[0].append(iele)

                    connect = [int(x)-1 for x in data[3+ntags:]]
                    self.addElement(connect)
        fid.close()

        # Print nnodes, nele and ngroups
        if self.ngroups == 1:
            logging.info(("Mesh read with {} nodes, {} elements and 1 group: ").format(
                        self.nnod, self.nele))
        else:
            logging.info(("Mesh read with {} nodes, {} elements and {} groups: ").format(
                self.nnod, self.nele, self.ngroups))

        # Print group names and number of elements
        for key in self.groupNames:
            group_name = self.groupNames[key]
            idx = self.groupNames.keys().index(key)
            group_nele = len(self.groups[idx])
            logging.info(("    {} with {} elements").format(group_name, group_nele))
        

    #-----------------------------------------------------------------------
    #   readXML
    #-----------------------------------------------------------------------

    def readXML(self, path=None, rank=3, doElemGroups=False):
        """ Input: self.path = self.path_to_file """
        if path is None:
            Tk().withdraw()
            self.path = filedialog.askopenfilename()
        else:
            self.path = path

        if doElemGroups is True:
            raise NotImplementedError(" readXML does not support doElemGroups!")
        else:
            self.groups = [[]]
            self.ngroups = 1
            self.groupNames[0] = 'Group 0'

        with open(self.path, 'r') as file:

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

                data = re.findall(r"[-+]?\d+ *\.\d+|[-+]?\d+", line)

                if len(data) > 0 and data[0].isdigit():

                    #-------------------------------------------------------
                    #   Nodes
                    #-------------------------------------------------------

                    if flag_n is True:
                        coord = [float(x) for x in data[1:rank+1]]
                        self.addNode(coord)

                    #-------------------------------------------------------
                    #   Elements
                    #-------------------------------------------------------

                    if flag_e is True:
                        connect = [int(x)-1 for x in data[1:]]
                        self.groups[0].append(int(data[0])-1)
                        self.addElement(connect)

        logging.info(("Mesh read with {} nodes, {} elements.").format(
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
            iele = connect + [connect[0]]
            coords = self.getCoords(iele)
            ax.plot(coords[:, 0], coords[:, 1], linewidth=0.5, color='k')

        return ax

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
            if idofs: # idofs is not empty
                x = self.getCoords(inod)
                u = np.array(disp[idofs])*scale
                deformed[inod, :] = u+x

        # Create figure
        if rank == 1 or rank == 2:
            fig = plt.figure(figsize=(6, 6))
            ax = fig.add_subplot(111)
        elif rank == 3:
            fig = plt.figure(figsize=(6, 6))
            ax = fig.add_subplot(111, projection='3d')

        # Plot deformed mesh
        for connect in self.connectivity:
            iele = connect + [connect[0]]
            coords = deformed[np.ix_(iele), :][0]
            ax.plot(coords[:, 0], coords[:, 1], linewidth=0.5, color='k')

        plt.show()

    #-----------------------------------------------------------------------
    #   updateGeometry
    #-----------------------------------------------------------------------

    def updateGeometry(self, disp):
        """ Input:  disp = displacement vector """

        for inod in range(len(self.coords)):
            idofs = self.getDofIndices(inod)
            if idofs:
                x = [a+b for a, b in zip(self.coords[inod], disp[idofs])]
                self.coords[inod] = x


#===========================================================================
#   Example
#===========================================================================


if __name__ == '__main__':

    mesh = Mesh()
    mesh.readGmsh("Examples/rve.msh")
    ax = mesh.plotMesh(rank=2)
    plt.show()
