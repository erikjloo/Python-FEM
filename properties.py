# Import Standard Libraries
import scipy as np
from enum import IntEnum


#===========================================================================
#   ElementType
#===========================================================================


class ElementType(IntEnum):
    # First order elements
    line2, tri3, quad4, tetra4, cube8, prism6, pyramid5 = range(1,8)
    # Second order elements
    line3, tri6, quad9, tetra10, cube27, prism18 = range(8,14)
    # Additional elements
    pyramid14, point1, quad8, cube20, prism15, pyramid13 = range(14,20)

    @classmethod
    def getTypeName(cls, etype):
        return cls(etype).name


#===========================================================================
#   Properties
#===========================================================================


class Properties(object):
    
    """ Properties

    Static Members:
        __type__ = "Input is not list or array!"
        __type_str__ = "Input is not str!"

    Instance Members:
        nele = number of elements
        materials = list of Materials
        types = list of property type names
        properties = array of property indices
            etype = properties[iele,0]
            imat = properties[iele,1]

    Public Methods:
        Material Methods:
            addMaterial(name, **kwargs)
            setMaterial(imat, name, **kwargs)
            eraseMaterial(name)

        Type Methods:
            addType(prop)
            addTypes(props)
            setType(jtype, prop)
            eraseType(prop)
            eraseTypes(props)
            ntyp = typeCount()
            prop = getTypeName(jtype)

        Property Methods:
            setProperty(iele, "mat", name)
            setProperty(iele, "etype", type)
            setProperty(iele, "other", value)
            [etype, mat, ...] = getProperties(iele,["etype","mat",...])

    Private Methods:
        __attachMaterial(iele, name)
        
    """
    # Static:
    __type__ = "Input is not list or array!"
    __type_str__ = "Input is not str!"

    class Material():
        def __init__(self, name, **kwargs):
            self.name = name
            self.__dict__.update(kwargs)
    
    # Public:
    def __init__(self, nele, properties=None):
        self.nele = nele
        self.materials = []
        self.types = ["etype","mat"]
        if properties is None:
            self.properties = np.zeros((nele,6))
            self.properties[:] = np.nan
        else:
            self.properties = properties

    #-------------------------------------------------------------------
    #   materials
    #-------------------------------------------------------------------

    def addMaterial(self, name, **kwargs):
        """ Input: name = string of material name, **kwargs """
        if isinstance(name, str):
            mat = self.Material(name, **kwargs)
            self.materials.append(mat)
        else:
            raise TypeError(self.__type_str__)

    def setMaterial(self, imat, name, **kwargs):
        """ Input: imat = material number, name = string of material name, **kwargs """
        if isinstance(name, str):
            mat = self.Material(name, **kwargs)
            self.materials[imat] = mat
        else:
            raise TypeError(self.__type_str__)
    
    def eraseMaterial(self, name):
        """ Input: name = string of material name to be erased """
        for imat, mat in enumerate(self.materials):
            if mat.name == name:
                del self.materials[imat]
    
    def eraseMaterials(self, names):
        """ Input: names = list of string of material names to be erased """
        for name in names:
            self.eraseMaterial(name)

    #-------------------------------------------------------------------
    #   types vector
    #-------------------------------------------------------------------

    def addType(self, prop):
        """ Input: prop = string of property type name """
        if isinstance(prop, str):
            self.types.append(prop)
        else:
            raise TypeError(self.__type_str__)
        # Check if dofspace has enough columns for all dof types
        c_reqd = len(self.types) - np.size(self.properties, 1)
        if c_reqd > 0:
            c_new = np.empty((self.nele, c_reqd))
            c_new[:] = np.nan
            self.properties = np.c_[self.properties, c_new]

    def addTypes(self, props):
        """ Input: props =  list of string of property type names """
        if isinstance(props, (list, np.ndarray)):
            for prop in props:
                self.addType(prop)
        else:
            raise TypeError(self.__type__)

    def setType(self, jtype, prop):
        """ Input: jtype = prop type index, prop = string of property type name """
        if isinstance(prop, str) and jtype > 1:
            self.types[jtype] = prop
        else:
            raise TypeError(self.__type_str__)

    def eraseType(self, prop):
        """ Input: prop = string of property type name to be erased """
        # Delete column from dofspace
        if prop != "etype" and prop != "mat":
            jtype = self.types.index(prop)
            self.properties = np.delete(self.properties, jtype, 1)
            del self.types[jtype]

    def eraseTypes(self, props):
        """ Input: props = list of strings of property type names to be erased """
        for prop in props:
            self.eraseType(prop)

    def typeCount(self):
        """ Output: number of property types """
        return len(self.types)

    def getTypeName(self, jtype):
        """ Input: jtype = property type index
            Output: string of property type name """
        return self.types[jtype]

    #-------------------------------------------------------------------
    #   properties
    #-------------------------------------------------------------------

    def setProperty(self, iele, prop, tag):
        """ Input:  iele = element index
                    prop = "etype", "mat", or user-specified property 
                    tag = etype number, material name, or user value """
        if prop == "etype":
            self.properties[iele, 0] = tag
        elif prop == "mat":
            self.__attachMaterial(iele, tag)
        else:
            jtype = self.types.index(prop)
            self.properties[iele, jtype] = tag

    def getProperties(self, iele, props):
        """ Input:  iele = element index
                    props = list of property type names
            Output: specified by props """
        output = []
        for prop in props:
            if prop == "mat":
                imat = int(self.properties[iele][1])
                output.append(self.materials[imat])
            else:
                jtype = self.types.index(prop)
                output.append(self.properties[iele][jtype])
        return output

    # Private:
    def __attachMaterial(self, iele, name):
        """ Input: iele = element index, name = material name """
        for imat, mat in enumerate(self.materials):
            if mat.name == name:
                self.properties[iele, 1] = imat
    

#===========================================================================
#   Example
#===========================================================================


if __name__ == "__main__":

    props = Properties(10)
    props.addMaterial("Matrix", E = 9000, v = 0.3)
    props.addMaterial("Matrix", E=9000, v=0.3)
    props.addMaterial("Matrix", E=9000, v=0.3)
    props.setMaterial(1,"Fibers", E = 45000, v = 0.2)
    props.addMaterial("Fibers2", E=9000, v=0.3)
    print("\n")
    [print(props.materials[i].name) for i in range(4)]

    props.eraseMaterial("Matrix")
    print("\n")
    [print(props.materials[i].name) for i in range(2)]

    print("\n")

    props.setProperty(0,"etype",2)
    props.setProperty(0,"mat","Fibers")
    # pylint: disable = unbalanced-tuple-unpacking
    [mat,etype] = props.getProperties(0,["mat","etype"])
    print(mat.name)
    print(etype)
