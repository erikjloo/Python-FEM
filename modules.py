#  Import Standard libraties
import scipy as np
from enum import IntEnum
from copy import deepcopy
from warnings import warn

#  Import Local Libraries
from abc import ABCMeta, abstractmethod
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
        for module in self.modules[:]:
            status = module.init(props, globdat)
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
        init(props, globdat)
        run(globdat)
        shutdown(globdat)
    """

    def init(self, props, globdat):
        props = props.getProps(self.name)
        globdat.makeMesh(props)
        return Status.DONE

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
        return Status.DONE

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
        try:
            props = props.getProps(self.name)
            self.type = props.get("solver.type")
        except KeyError:
            warn("No solver specified ")
            self.type = "solve"
        return Status.OK

    def run(self, globdat):
        """ Runs a linear analysis """

        # ADVANCE
        print("Advancing to next load step")
        # globdat.model.advance()

        # GET_MATRIX_0
        self.hbw = globdat.model.get_Matrix_0(
            globdat.mbuild, globdat.fint, globdat.disp, globdat.mesh)
        print("Stiffness matrix: hbw = {} out of {}".format(
            self.hbw, globdat.mesh.dofCount()))

        # GET_EXT_VECTOR
        globdat.model.get_Ext_Vector(globdat.fext, globdat.mesh)
        globdat.fext += globdat.load.getLoads()
        print("External force vector: {} elements".format(len(globdat.fext)))

        # GET_CONSTRAINTS
        globdat.model.get_Constraints(globdat.cons, globdat.mesh)
        globdat.disp = globdat.cons.getDisps()
        sdof = globdat.cons.get_sdof()
        print("Constraints: {} prescribed values".format(len(sdof)))

        # INITIAL RESIDUAL
        self.solver = Solver(self.type, globdat.cons)
        K = globdat.mbuild.getDenseMatrix()
        r = globdat.fext - globdat.fint - np.array(K).dot(globdat.disp)

        self.solver.solve(K, globdat.disp, r, self.hbw)

        # GET_INT_VECTOR
        globdat.model.get_Int_Vector(globdat.fint, globdat.disp, globdat.mesh)

        return Status.DONE

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
        """ Sets up solver """
        
        try:
            props = props.getProps(self.name)
            self.nrkey = props.get("type")
            self.niter = props.get("niter")
            self.tol = props.get("tol")
            self.type = props.get("solver.type")
        except KeyError:
            warn("No solver specified ")
            self.nrkey = "full"
            self.niter = 20
            self.tol = 1e-3
            self.type = "solve"

        # GET_MATRIX_0
        self.hbw = globdat.model.get_Matrix_0(
            globdat.mbuild, globdat.fint, globdat.disp, globdat.mesh)
        print("Stiffness matrix: hbw = {} out of {}".format(
            self.hbw, globdat.mesh.dofCount()))

    def run(self, globdat):
        """ Runs a non-linear analysis """
        
        # ADVANCE
        print("Advancing to next load step")
        # globdat.model.advance()

        # GET_EXT_VECTOR
        globdat.model.get_Ext_Vector(globdat.fext, globdat.mesh)
        globdat.fext += globdat.load.getLoads()
        print("External force vector: {} elements".format(len(globdat.fext)))

        # GET_CONSTRAINTS
        globdat.model.get_Constraints(globdat.cons, globdat.mesh)
        Du = globdat.cons.getDisps()
        du = deepcopy(Du)
        sdof = globdat.cons.get_sdof()
        globdat.disp[sdof] += Du[sdof]
        print("Constraints: {} prescribed values".format(len(sdof)))

        # INITIAL RESIDUAL
        self.solver = Solver(self.type, globdat.cons)
        K = globdat.mbuild.getDenseMatrix()
        r = globdat.fext - globdat.fint - np.array(K).dot(Du)
        fdof = globdat.cons.get_fdof()
        # du = np.zeros(globdat.ndof)

        for iter in range(self.niter):
            
            # Update displacement vector
            self.solver.solve(K, du, r, self.hbw)
            # Du[fdof] += du[fdof]
            globdat.disp[fdof] += du[fdof]

            # Find interal force vector
            if self.nrkey == "full":
                globdat.model.get_Matrix_0(
                    globdat.mbuild, globdat.fint, globdat.disp, globdat.mesh)
            elif self.nrkey == "mod" or self.nrkey == "LE":
                globdat.model.get_Int_Vector(
                    globdat.fint, du, globdat.mesh)

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
