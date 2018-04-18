
# Import Standard Libraries
import scipy as np
from abc import ABC, abstractmethod


#===========================================================================
#   Material
#===========================================================================


class Material(ABC):

	def __init__(self, props):

    	for name, val in props:
      		setattr(self, name, val)

    self.initHistory = {}
    self.current = []
    self.iIter = -1

  	def setIter(self, iIter):
		self.iIter = iIter

  	def setHistoryParameter(self, name, val):

    	if self.iIter == -1:
      		self.initHistory[name] = val
      		return

    	if len(self.current) == self.iIter:
      		self.current.append(self.initHistory.copy())

    	self.current[self.iIter][name] = val
			return

		def getHistoryParameter(self, name):

    		if len(self.history) == 0:
      			return self.initHistory[name]
    		else:
      			return self.history[self.iIter][name]

  		def commitHistory(self):

    		self.history = []

    		for h in self.current:
      			self.history.append(h)


#===========================================================================
#   Mesh
#===========================================================================


class PlaneStrain( BaseMaterial ):

  def __init__ ( self, props ):

    #Call the BaseMaterial constructor
    BaseMaterial.__init__( self, props )

    #Create the hookean matrix
    self.H = np.zeros( (3,3) )

    self.H[0,0] = self.E*(1.-self.nu)/((1+self.nu)*(1.-2.*self.nu));
    self.H[0,1] = self.H[0,0]*self.nu/(1-self.nu);
    self.H[1,0] = self.H[0,1];
    self.H[1,1] = self.H[0,0];
    self.H[2,2] = self.H[0,0]*0.5*(1.-2.*self.nu)/(1.-self.nu);

  def getStress( self, deformation ):

    sigma = dot( self.H, deformation.strain )

    return sigma, self.H

  def getTangent( self ):
  
    return self.H
