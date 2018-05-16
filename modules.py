#  Import Standard libraties
import scipy as np

#  Import Local Libraries
from abc import ABCMeta, abstractmethod
from solvers import Solver


#===========================================================================
#   Module
#===========================================================================


class Module(metaclass=ABCMeta):
    """ Abstract Module Class 
    
    Pure Virtual Methods:
        init(props, globdat)
        run(globdat)
        shutdown(globdat)
    """

    def __init__(self, name=None):
        self.name = name

    @abstractmethod
    def init(self, props, globdat): pass

    @abstractmethod
    def run(self, globdat): pass

    @abstractmethod
    def shutdown(self, globdat): pass


#===========================================================================
#   ChainModule
#===========================================================================


class ChainModule(Module):
    """ Groups several modules in order of execution
        
    Public Methods:
        pushBack(module)
        init(props, globdat)
        run(globdat)
        shutdown(globdat)
    """

    def __init__(self, name=None):
        self.name = name
        self.modules = []

    def pushBack(self, module):
        self.modules.append(module)

    def init(self, props, globdat):
        for module in self.modules:
            module.init(props, globdat)

    def run(self, globdat):
        for module in self.modules:
            module.run(globdat)

    def shutdown(self, globdat):
        for module in self.modules:
            module.shutdown(globdat)


#===========================================================================
#   InputModule
#===========================================================================


class InputModule(Module):
    """ Reads props and initializes the mesh
        
    Public Methods:
        init(props, globdat)
        run(globdat)
        shutdown(globdat)
    """

    def init(self, props, globdat):
        props = props.getProps(self.name)
        globdat.makeMesh(props)

    def run(self, globdat): pass

    def shutdown(self, globdat): pass


#===========================================================================
#   InitModule
#===========================================================================


class InitModule(Module):
    """ Initializes the model, loads, constraints, matrix builder & vectors
            
    Public Methods:
        init(props, globdat)
        run(globdat)
        shutdown(globdat)
    """

    def init(self, props, globdat):
        globdat.makeModel(props)
        globdat.makeLoadTable(props)
        globdat.makeConstraints(props)
        globdat.makeMatrixBuilder()
        globdat.makeVectors()

    def run(self, globdat): pass

    def shutdown(self, globdat): pass


#===========================================================================
#   LinSolveModule
#===========================================================================


class LinSolveModule(Module):
    """ Runs a linear analysis
            
    Public Methods:
        init(props, globdat)
        run(globdat)
        shutdown(globdat)
    """

    def init(self, props, globdat):
        """ Computes K, f_ext and cons """

        print("    Assembling stiffness matrix")
        self.hbw = globdat.model.get_Matrix_0(
            globdat.mbuild, globdat.fint, globdat.disp, globdat.mesh)
        print("        The half-band-width is", self.hbw)
        
        print("    Assembling external force vector")
        globdat.fext = globdat.load.getLoads()
        globdat.model.get_Ext_Vector(globdat.fext, globdat.mesh)

        print("    Updating constraints")
        globdat.model.get_Constraints(globdat.cons, globdat.mesh)

    def run(self, globdat):
        """ Runs a linear analysis """

        print("    Running linear analysis")
        solver = Solver("numpy", globdat.cons)
        K = globdat.mbuild.getDenseMatrix()
        solver.solve(K, globdat.disp, globdat.fext, self.hbw)
        
    def shutdown(self, globdat):
        pass


#===========================================================================
#   NonlinModule
#===========================================================================


class NonlinModule(Module):
    """ Runs a nonlinear analysis (NR)
            
    Public Methods:
        init(props, globdat)
        run(globdat)
        shutdown(globdat)
    """

    def init(self, props, globdat):
        """ Computes K, f_ext and cons """

        print("    Assembling stiffness matrix")
        self.hbw = globdat.model.get_Matrix_0(
            globdat.mbuild, globdat.fint, globdat.disp, globdat.mesh)
        print("        The half-band-width is", self.hbw)

        print("    Assembling external force vector")
        globdat.fext = globdat.load.getLoads()
        globdat.model.get_Ext_Vector(globdat.fext, globdat.mesh)

        print("    Updating constraints")
        globdat.model.get_Constraints(globdat.cons, globdat.mesh)

    def run(self, globdat):
        """ Runs a linear analysis """

        print("    Running linear analysis")
        solver = Solver("numpy", globdat.cons)
        K = globdat.mbuild.getDenseMatrix()
        solver.solve(K, globdat.disp, globdat.fext, self.hbw)

    def shutdown(self, globdat):
        pass
