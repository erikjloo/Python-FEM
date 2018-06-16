#  Import Standard libraties
import scipy as np
from enum import IntEnum
from copy import deepcopy
from warnings import warn

#  Import Local Libraries
from abc import ABCMeta, abstractmethod
from properties import Properties
from algebra import norm
from solvers import Solver


#===========================================================================
#   Status
#===========================================================================


class Status(IntEnum):
    OK = 0
    DONE = 1
    EXIT = 2


#===========================================================================
#   Module
#===========================================================================


class Module(metaclass=ABCMeta):
    """ Abstract Module Class 
    
    Pure Virtual Methods:
        init(props, conf, globdat)
        run(globdat)
        shutdown(globdat)
    """

    def __init__(self, name=None):
        self.name = name

    @abstractmethod
    def init(self, props, conf, globdat): pass

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
        init(props, conf, globdat)
        run(globdat)
        shutdown(globdat)
    """

    def __init__(self, name=None):
        self.name = name
        self.modules = []

    def pushBack(self, module):
        self.modules.append(module)

    def init(self, props, conf, globdat):
        for module in self.modules[:]:
            status = module.init(props, conf, globdat)
            if status == Status.DONE:
                self.modules.remove(module)
            elif status == Status.EXIT:
                break

    def run(self, globdat):
        for module in self.modules[:]:
            status = module.run(globdat)
            if status == Status.DONE:
                self.modules.remove(module)
            elif status == Status.EXIT:
                break
        return status

    def shutdown(self, globdat):
        for module in self.modules:
            module.shutdown(globdat)


#===========================================================================
#   InputModule
#===========================================================================


class InputModule(Module):
    """ Reads props and initializes the mesh
        
    Public Methods:
        init(props, conf, globdat)
        run(globdat)
        shutdown(globdat)
    """

    def init(self, props, conf, globdat):
        myProps = props.getProps(self.name)
        myConf = conf.makeProps(self.name)
        globdat.makeMesh(myProps,myConf)
        return Status.DONE

    def run(self, globdat): pass

    def shutdown(self, globdat): pass


#===========================================================================
#   InitModule
#===========================================================================


class InitModule(Module):
    """ Initializes the model, loads, constraints, matrix builder & vectors
            
    Public Methods:
        init(props, conf, globdat)
        run(globdat)
        shutdown(globdat)
    """

    def init(self, props, conf, globdat):
        globdat.makeModel(props, conf)
        globdat.makeLoadTable(props)
        globdat.makeConstraints(props)
        globdat.makeMatrixBuilder()
        globdat.makeVectors()
        return Status.DONE

    def run(self, globdat): pass

    def shutdown(self, globdat): pass


#===========================================================================
#   LinSolveModule
#===========================================================================


class LinSolveModule(Module):
    """ Runs a linear analysis
            
    Public Methods:
        init(props, conf, globdat)
        run(globdat)
        shutdown(globdat)
    """

    def init(self, props, conf, globdat):
        """ Computes K, f_ext and cons """

        myProps = props.getProps(self.name)
        myConf = conf.makeProps(self.name)
        self.type = myProps.get("solver.type","solve")
        myConf.set("solver.type",self.type)
        return Status.OK

    def run(self, globdat):
        """ Runs a linear analysis """

        # ADVANCE
        print("Advancing to next load step")
        # globdat.model.advance()

        # GET_MATRIX_0
        globdat.model.takeAction("GET_MATRIX_0",globdat)
        print("Stiffness matrix: hbw = {} out of {}".format(
            globdat.mbuild.hbw, globdat.mesh.dofCount()))

        # GET_EXT_VECTOR
        globdat.model.takeAction("GET_EXT_VECTOR",globdat)
        globdat.fext += globdat.load.getLoads()
        print("External force vector: {} elements".format(len(globdat.fext)))

        # GET_CONSTRAINTS
        globdat.model.takeAction("GET_CONSTRAINTS",globdat)
        globdat.disp = globdat.cons.getDisps()
        sdof = globdat.cons.get_sdof()
        print("Constraints: {} prescribed values".format(len(sdof)))

        # INITIAL RESIDUAL
        self.solver = Solver(self.type, globdat.cons)
        K = globdat.mbuild.getDenseMatrix()
        r = globdat.fext - globdat.fint - np.array(K).dot(globdat.disp)

        self.solver.solve(K, globdat.disp, r, globdat.mbuild.hbw)

        # GET_INT_VECTOR
        globdat.model.takeAction("GET_INT_VECTOR", globdat)

        return Status.DONE

    def shutdown(self, globdat):
        pass


#===========================================================================
#   NonlinModule
#===========================================================================


class NonlinModule(Module):
    """ Runs a nonlinear analysis (NR)
            
    Public Methods:
        init(props, conf, globdat)
        run(globdat)
        shutdown(globdat)
    """

    def init(self, props, conf, globdat):
        """ Sets up solver """
        
        myProps = props.getProps(self.name)
        myConf = conf.makeProps(self.name)

        self.nrkey = myProps.get("type","full")
        self.niter = myProps.get("niter",2)
        self.tol = myProps.get("tol",1e-3)
        self.type = myProps.get("solver.type","solve")

        myConf.set("type",self.nrkey)
        myConf.set("niter",self.niter)
        myConf.set("tol",self.tol)
        myConf.set("solver.type",self.type)

    def run(self, globdat):
        """ Runs a non-linear analysis """
        
        # ADVANCE
        print("Advancing to next load step")
        # globdat.model.advance()
        
        # GET_MATRIX_0
        globdat.model.takeAction("GET_MATRIX_0", globdat)
        print("Stiffness matrix: hbw = {} out of {}".format(
            globdat.mbuild.hbw, globdat.mesh.dofCount()))

        # GET_EXT_VECTOR
        globdat.fext = np.zeros(globdat.ndof)
        globdat.model.takeAction("GET_EXT_VECTOR", globdat)
        globdat.fext += globdat.load.getLoads()
        print("External force vector: {} elements".format(len(globdat.fext)))

        # GET_CONSTRAINTS
        globdat.model.takeAction("GET_CONSTRAINTS", globdat)
        du = globdat.cons.getDisps()
        sdof = globdat.cons.get_sdof()
        disp_old = deepcopy(globdat.disp)
        globdat.disp[sdof] += du[sdof]
        print("Constraints: {} prescribed values".format(len(sdof)))

        # INITIAL RESIDUAL
        self.solver = Solver(self.type, globdat.cons)
        K = globdat.mbuild.getDenseMatrix()
        r = globdat.fext - globdat.fint - np.array(K).dot(du)
        fdof = globdat.cons.get_fdof()

        for iter in range(self.niter):
            
            # Update displacement vector
            self.solver.solve(K, du, r, globdat.mbuild.hbw)
            globdat.disp[fdof] += du[fdof]

            # Find interal force vector
            globdat.fint = np.zeros(globdat.ndof)
            
            if self.nrkey == "full":
                globdat.model.takeAction("GET_MATRIX_0",globdat)
            elif self.nrkey == "mod" or self.nrkey == "LE":
                globdat.model.takeAction("GET_INT_VECTOR",globdat)

            # Find out-of-balance force vector
            r = globdat.fext - globdat.fint
            nrm = norm(r[fdof])
            print("    Iteration {}: norm = {:.2f} ".format(iter,nrm))
            
            # Check convergence
            if iter == 0:
                nrm1 = deepcopy(nrm)
            if nrm < self.tol*nrm1:
                print("    Converged in {} iterations".format(iter+1))
                break
        
        # globdat.mesh.updateGeometry(Du)
        # Commit
        
    def shutdown(self, globdat):
        pass


#===========================================================================
#   NonlinModule
#===========================================================================


class OutputModule(Module):
    """ Saves specified output
            
    Public Methods:
        init(props, conf, globdat)
        run(globdat)
        shutdown(globdat)
    """

    def init(self, props, conf, globdat):
        myProps = props.getProps(self.name)
        myConf = props.makeProps(self.name)
        return Status.DONE

    def run(self, globdat): pass

    def shutdown(self, globdat): pass
