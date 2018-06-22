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
        Status = init(conf, props, globdat)
        Status = run(globdat)
        shutdown(globdat)
    """

    def __init__(self, name=None):
        self.name = name

    def __del__(self):
        print("Cleaning {} module".format(self.name))

    @abstractmethod
    def init(self, conf, props, globdat):
        raise NotImplementedError()

    @abstractmethod
    def run(self, globdat):
        raise NotImplementedError()

    @abstractmethod
    def shutdown(self, globdat):
        raise NotImplementedError()


#===========================================================================
#   ChainModule
#===========================================================================


class ChainModule(Module):
    """ Groups several modules in order of execution

    Public Methods:
        pushBack(module)
        Status = init(conf, props, globdat)
        Status = run(globdat)
        shutdown(globdat)
    """

    def __init__(self, name="chain"):
        self.name = name
        self.modules = []

    def pushBack(self, module):
        self.modules.append(module)

    def init(self, conf, props, globdat):
        for module in self.modules[:]:
            status = module.init(conf, props, globdat)
            if status == Status.DONE:
                self.modules.remove(module)
            elif status == Status.EXIT:
                return Status.EXIT
        return Status.OK

    def run(self, globdat):
        for module in self.modules[:]:
            status = module.run(globdat)
            if status == Status.DONE:
                self.modules.remove(module)
            elif status == Status.EXIT:
                return status.EXIT
        return status.OK

    def shutdown(self, globdat):
        for module in self.modules:
            module.shutdown(globdat)


#===========================================================================
#   ControlModule
#===========================================================================


class ControlModule(Module):
    """ Controls the number of steps

    Public Methods:
        Status = init(conf, props, globdat)
        Status = run(globdat)
        shutdown(globdat)
    """

    __key__ = "Control: nsteps not specified!"
    
    def __init__(self, name="control"):
        self.name = name

    def init(self, conf, props, globdat):
        myProps = props.getProps(self.name)
        myConf = conf.makeProps(self.name)
        
        self.nsteps = myProps.get("nsteps")
        myConf.set("nsteps",self.nsteps)

        if self.nsteps is None:
            raise KeyError(self.__key__)

        return Status.OK

    def run(self, globdat):
        if globdat.i < self.nsteps:
            return Status.OK
        elif globdat.i == self.nsteps:
            return Status.EXIT

    def shutdown(self, globdat):
        pass


#===========================================================================
#   InputModule
#===========================================================================


class InputModule(Module):
    """ Reads props and initializes the mesh
        
    Public Methods:
        Status = init(conf, props, globdat)
        Status = run(globdat)
        shutdown(globdat)
    """

    def __init__(self, name="input"):
        self.name = name

    def init(self, conf, props, globdat):
        myProps = props.getProps(self.name)
        myConf = conf.makeProps(self.name)
        globdat.makeMesh(myProps, myConf)
        return Status.DONE

    def run(self, globdat): 
        pass

    def shutdown(self, globdat): 
        pass


#===========================================================================
#   InitModule
#===========================================================================


class InitModule(Module):
    """ Initializes the model, loads, constraints, matrix builder & vectors
            
    Public Methods:
        Status = init(conf, props, globdat)
        Status = run(globdat)
        shutdown(globdat)
    """

    def __init__(self, name="init"):
        self.name = name

    def init(self, conf, props, globdat):
        globdat.makeModel(props, conf)
        globdat.makeLoadTable(props)
        globdat.makeConstraints(props)
        globdat.makeMatrixBuilder()
        globdat.makeVectors()
        return Status.DONE

    def run(self, globdat):
        pass

    def shutdown(self, globdat): 
        pass


#===========================================================================
#   LinSolveModule
#===========================================================================


class LinSolveModule(Module):
    """ Runs a linear analysis
            
    Public Methods:
        Status = init(conf, props, globdat)
        Status = run(globdat)
        shutdown(globdat)
    """

    def __init__(self, name="linsolve"):
        self.name = name

    def init(self, conf, props, globdat):
        myProps = props.getProps(self.name)
        myConf = conf.makeProps(self.name)

        self.type = myProps.get("solver.type", "solve")
        myConf.set("solver.type", self.type)

        return Status.OK

    def run(self, globdat):

        # ADVANCE
        print("Advancing to the next load step")
        globdat.i += 1
        globdat.model.takeAction("ADVANCE", globdat)

        # GET_MATRIX_0
        globdat.model.takeAction("GET_MATRIX_0", globdat)
        print("Stiffness matrix: {} dofs".format(globdat.ndof))

        # GET_EXT_VECTOR
        globdat.model.takeAction("GET_EXT_VECTOR", globdat)
        globdat.fext += globdat.load.getLoads()
        print("External force vector: {} dofs".format(len(globdat.fext)))

        # GET_CONSTRAINTS
        globdat.model.takeAction("GET_CONSTRAINTS", globdat)
        globdat.disp = globdat.cons.getDisps()
        sdof = globdat.cons.get_sdof()
        print("Constraints: {} prescribed values".format(len(sdof)))

        # INITIAL RESIDUAL
        self.solver = Solver(self.type, globdat.cons)
        K = globdat.mbuild.getDenseMatrix()
        r = globdat.fext - globdat.fint - np.array(K).dot(globdat.disp)

        self.solver.solve(K, globdat.disp, r, globdat.mbuild.hbw)

        # GET_INT_VECTOR
        globdat.fint = np.zeros(globdat.ndof)
        globdat.model.takeAction("GET_INT_VECTOR", globdat)

        return Status.EXIT

    def shutdown(self, globdat): 
        pass


#===========================================================================
#   NonlinModule
#===========================================================================


class NonlinModule(Module):
    """ Runs a nonlinear analysis (NR)
            
    Public Methods:
        Status = init(conf, props, globdat)
        Status = run(globdat)
        shutdown(globdat)
    """

    def __init__(self, name="nonlin"):
        self.name = name

    def init(self, conf, props, globdat):
        myProps = props.getProps(self.name)
        myConf = conf.makeProps(self.name)

        self.nrkey = myProps.get("type", "full")
        self.niter = myProps.get("niter",20)
        self.tiny = myProps.get("tiny",1e-10)
        self.tol = myProps.get("tol", 1e-3)
        self.type = myProps.get("solver.type", "solve")

        myConf.set("type", self.nrkey)
        myConf.set("niter", self.niter)
        myConf.set("tiny", self.tiny)
        myConf.set("tol", self.tol)
        myConf.set("solver.type", self.type)

        return Status.OK

    def run(self, globdat):

        # ADVANCE
        print("Advancing to the next load step")
        globdat.i += 1
        globdat.model.takeAction("ADVANCE", globdat)

        # GET_MATRIX_0
        globdat.model.takeAction("GET_MATRIX_0", globdat)
        print("Stiffness matrix: {} dofs".format(globdat.ndof))

        # GET_EXT_VECTOR
        globdat.fext = np.zeros(globdat.ndof)
        globdat.model.takeAction("GET_EXT_VECTOR", globdat)
        print("External force vector: {} dofs".format(len(globdat.fext)))

        # GET_CONSTRAINTS
        globdat.model.takeAction("GET_CONSTRAINTS", globdat)
        du = globdat.cons.getDisps()
        sdof = globdat.cons.get_sdof()
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
                globdat.model.takeAction("GET_MATRIX_0", globdat)
            elif self.nrkey == "mod" or self.nrkey == "LE":
                globdat.model.takeAction("GET_INT_VECTOR", globdat)
            else:
                raise ValueError("{} not implemented !".format(self.nrkey))

            # Find out-of-balance force vector
            r = globdat.fext - globdat.fint
            nrm = norm(r[fdof])
            print("    Iteration {}: norm = {:.10f} ".format(iter, nrm))
            
            # Check convergence in first iteration
            if iter == 0 and nrm <= self.tiny:
                print("    Converged in {} iterations".format(iter+1))
                return Status.OK
            elif iter == 0 and nrm > self.tiny:
                nrm1 = deepcopy(nrm)
            
            # Check convergence in later iterations
            if nrm < self.tol*nrm1:
                print("    Converged in {} iterations".format(iter+1))
                return Status.OK

    def shutdown(self, globdat): 
        pass


#===========================================================================
#   SampleModule
#===========================================================================


class SampleModule(Module):
    """ Samples force-displacement points
            
    Public Methods:
        Status = init(conf, props, globdat)
        Status = run(globdat)
        shutdown(globdat)
    """

    __key__ = "Sample: file or dofs not specified!"

    def __init__(self, name="sample"):
        self.name = name

    def init(self, conf, props, globdat):
        myProps = props.getProps(self.name)
        myConf = conf.makeProps(self.name)

        self.path = myProps.get("file")
        self.dofs = myProps.get("dofs")

        if self.path is None or self.dofs is None:
            raise KeyError(self.__key__)
        
        myConf.set("file", self.path)
        myConf.set("dofs", self.dofs)
        open(self.path, "w+").close()

        self.run(globdat)
        return Status.OK

    def run(self, globdat):
        i = globdat.i
        f = globdat.fint[self.dofs]
        u = globdat.disp[self.dofs]
        txt = "{}".format(i)

        if len(self.dofs) > 1:
            for u_f in zip(u, f):
                txt += " {} {} ".format(u_f[0], u_f[1])
            txt += " \n"
        elif len(self.dofs) == 1:
            txt += " {} {} ".format(u, f)

        with open(self.path, 'a') as f:
            f.write(txt)
        return Status.OK

    def shutdown(self, globdat): 
        pass 


#===========================================================================
#   Execute
#===========================================================================


def Execute(module, conf, props, globdat):

    # Initialize all modules
    print("==== Initializing modules ================")
    status = module.init(conf, props, globdat)
    conf.print()
    globdat.model.takeAction("PLOT_MESH", globdat)
    globdat.mesh.printDofSpace(15)
    # Run modules until Status.EXIT is issued
    while status != Status.EXIT:
        print("==== Running modules =====================")
        status = module.run(globdat)

    # Shutdown remaining modules
    print("==== Shutting down modules ===============")
    module.shutdown(globdat)
