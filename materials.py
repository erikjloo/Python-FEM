# Import Standard Libraries
import scipy as np
from abc import ABCMeta, abstractmethod


#===========================================================================
#   MaterialFactory
#===========================================================================


def MaterialFactory(props):

    props = props.getProps("material")
    type = props.get("type")

    if type == "Hooke":
        print("Creating Hooke material")
        return Hooke(props)
    if type == "Melro":
        print("Creating Melro material")
        return PlaneStrain(props)
    elif type == "PlaneStrain":
        print("Creating PlaneStress material")
        return PlaneStrain(props)
    elif type == "PlaneStress":
        print("Creating PlaneStress material")
        return PlaneStress(props)


#===========================================================================
#   Material
#===========================================================================


class Material(metaclass=ABCMeta):
    """ Material """
    
    @abstractmethod
    def __init__(self, props):
        pass

    @abstractmethod
    def getStress(self, strain):
        pass
    
    @abstractmethod
    def getTangent(self):
        pass


#===========================================================================
#   Hooke
#===========================================================================


class Hooke(Material):

    def __init__(self, props):

        E = props.get("young")
        nu = props.get("poisson")

        # Calculate the Lamé parameters
        self.la = nu*E/((1+nu)*(1-2*nu))
        self.mu = E/(2*(1+nu))

        # Create the hookean matrix
        self.H = np.zeros((6, 6))
        self.H[np.ix_([0, 1, 2], [0, 1, 2])] = self.la
        self.H[0, 0] = self.H[1, 1] = self.H[2, 2] = self.la + 2*self.mu
        self.H[3, 3] = self.H[4, 4] = self.H[5, 5] = self.mu

    def getStress(self, strain):
        stress = self.H.dot(strain)
        return [stress, self.H]

    def getTangent(self):
        return self.H


#===========================================================================
#   PlaneStrain
#===========================================================================


class PlaneStrain(Material):

    def __init__(self, props):

        E = props.get("young")
        nu = props.get("poisson")

        # Calculate the Lamé parameters
        self.la = nu*E/((1+nu)*(1-2*nu))
        self.mu = E/(2*(1+nu))

        #Create the hookean matrix
        self.H = np.zeros((3, 3))
        self.H[0, 0] = self.H[1, 1] = self.la+2*self.mu
        self.H[0, 1] = self.H[1, 0] = self.la
        self.H[2, 2] = self.mu

    def getStress(self, strain):
        stress = self.H.dot(strain)
        return [stress, self.H]

    def getTangent(self):
        return self.H


#===========================================================================
#   PlaneStress
#===========================================================================


class PlaneStress(Material):

    def __init__(self, props):

        E = props.get("young")
        nu = props.get("poisson")

        # Calculate the Lamé parameters
        self.la = E*nu/(1-nu**2)
        self.mu = E/(2*(1+nu))

        #Create the hookean matrix
        self.H = np.zeros((3, 3))
        self.H[0, 0] = self.H[1, 1] = self.la+2*self.mu
        self.H[0, 1] = self.H[1, 0] = self.la
        self.H[2, 2] = self.mu

    def getStress(self, strain):
        stress = self.H.dot(strain)
        return [stress, self.H]

    def getTangent(self):
        return self.H


#===========================================================================
#   Melro
#===========================================================================


class Melro(Material):

    def __init__(self, props):

        E = props.get("young")
        nu = props.get("poisson")

        # Calculate the Lamé parameters
        self.la = nu*E/((1+nu)*(1-2*nu))
        self.mu = E/(2*(1+nu))

        # Create the hookean matrix
        self.H = np.zeros((6, 6))
        self.H[np.ix_([0, 1, 2], [0, 1, 2])] = self.la
        self.H[0, 0] = self.H[1, 1] = self.H[2, 2] = self.la + 2*self.mu
        self.H[3, 3] = self.H[4, 4] = self.H[5, 5] = self.mu

    def getStress(self, strain):
        stress = self.H.dot(strain)
        return [stress, self.H]

    def getTangent(self):
        return self.H
