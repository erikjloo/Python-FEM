#  Import Standard libraties
import scipy as np
from enum import IntEnum
from copy import deepcopy

#  Import Local Libraries
from abc import ABCMeta, abstractmethod
from mesh import Mesh
from models import Model, Action
from algebra import MatrixBuilder, norm
from loadTable import LoadTable
from constraints import Constraints
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

    def __init__(self, name="control"):
        self.name = name

    def init(self, conf, props, globdat):
        myProps = props.getProps(self.name)
        myConf = conf.makeProps(self.name)
        
        self.nsteps = myProps.get("nsteps")
        myConf.set("nsteps",self.nsteps)

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
    """ Reads props and initializes the mesh, loads and constraints
        
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

        modules = myProps.get("modules")
        myConf.set("modules", modules)

        for module in modules:
            type = myProps.get("{}.type".format(module))
            if type == "Gmsh" or type == "XML":
                mesh = Mesh(myConf, myProps)
                globdat.set("mesh",mesh)
            elif type == "Loads":
                load = LoadTable(module, myConf, myProps)
                globdat.set(module, load)
            elif type == "Constraints":
                cons = Constraints(module, myConf, myProps)
                globdat.set(module, cons)
            else:
                raise KeyError("Unknown input type: {}".format(type))
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
        globdat.model = Model.modelFactory("model", conf, props, globdat)
        ndof = globdat.set("ndof",globdat.get("mesh").dofCount())
        globdat.set("mbuild", MatrixBuilder(ndof))
        globdat.set("fext", np.zeros(ndof))
        globdat.set("fint", np.zeros(ndof))
        globdat.set("solu", np.zeros(ndof))
        globdat.set("loadScale",1.0)
        globdat.i = 0
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
        ndof = globdat.get("ndof")
        globdat.set("fext", np.zeros(ndof))
        globdat.set("loadScale", globdat.i)

        globdat.model.takeAction(Action.ADVANCE, globdat)
        globdat.model.takeAction(Action.GET_MATRIX_0, globdat)
        globdat.model.takeAction(Action.GET_EXT_VECTOR, globdat)
        globdat.model.takeAction(Action.GET_CONSTRAINTS, globdat)

        # INITIAL RESIDUAL
        mbuild = globdat.get("mbuild")
        fext = globdat.get("fext")
        cons = globdat.get("cons")
        disp = globdat.get("solu")
        ndof = globdat.get("ndof")

        old_disp = deepcopy(disp)
        cons.updateSolution(disp)
        Du = globdat.set("Du",disp - old_disp)

        K = mbuild.getDenseMatrix()
        r = fext - np.array(K).dot(Du)

        solver = Solver(self.type, cons)
        solver.solve(K, disp, r, mbuild.hbw)

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

        if self.nrkey not in ["NR","MNR","full","mod"]:
            raise ValueError("{} not implemented !".format(self.nrkey))
        if self.nrkey == "NR" or self.nrkey == "full":
            self.action = Action.GET_MATRIX_0 
        elif self.nrkey == "MNR" or self.nrkey == "mod":
            self.action = Action.GET_INT_VECTOR

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
        ndof = globdat.get("ndof")
        globdat.set("du", np.zeros(ndof))
        globdat.set("Du", np.zeros(ndof))
        globdat.set("fext", np.zeros(ndof))
        globdat.set("fint", np.zeros(ndof))
        globdat.set("loadScale", globdat.i)

        globdat.model.takeAction(Action.ADVANCE, globdat)
        globdat.model.takeAction(Action.GET_MATRIX_0, globdat)
        globdat.model.takeAction(Action.GET_EXT_VECTOR, globdat)
        globdat.model.takeAction(Action.GET_CONSTRAINTS, globdat)

        mbuild = globdat.get("mbuild")
        fext = globdat.get("fext")
        fint = globdat.get("fint")
        cons = globdat.get("cons")
        disp = globdat.get("solu")
        
        old_disp = deepcopy(disp)
        cons.updateSolution(disp)
        du = globdat.set("du", disp-old_disp)
        Du = globdat.set("Du", deepcopy(du))

        K = mbuild.getDenseMatrix()
        r = fext - fint - np.array(K).dot(Du)

        solver = Solver(self.type, cons)
        fdof = cons.getFdof()
        nrm1 = 0.0

        for iter in range(self.niter):

            # Update displacement vector
            solver.solve(K, du, r, mbuild.hbw)
            disp[fdof] += du[fdof]
            Du[fdof] += du[fdof]

            # Find interal force vector
            globdat.set("fint",np.zeros(ndof))
            globdat.model.takeAction(self.action, globdat)

            # Find out-of-balance force vector
            r = fext - globdat.get("fint")
            nrm = norm(r[fdof])
            print("    Iteration {}: norm = {:.10f} ".format(iter, nrm))
            
            # Check convergence
            if (iter == 0 and nrm <= self.tiny) or (iter > 0 and nrm < self.tol*nrm1):
                print("    Converged in {} iterations".format(iter+1))
                globdat.model.takeAction(Action.COMMIT, globdat)
                return Status.OK
            elif iter == 0 and nrm > self.tiny:
                nrm1 = deepcopy(nrm)
        
        return Status.EXIT

    def shutdown(self, globdat): 
        pass


#===========================================================================
#   ArclenModule
#===========================================================================


class ArclenModule(Module):
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
        self.niter = myProps.get("niter", 20)
        self.tiny = myProps.get("tiny", 1e-10)
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
        ndof = globdat.get("ndof")
        globdat.set("fext", np.zeros(ndof))
        globdat.set("loadScale", globdat.i)
        globdat.model.takeAction(Action.ADVANCE, globdat)
        globdat.model.takeAction(Action.GET_MATRIX_0, globdat)
        globdat.model.takeAction(Action.GET_EXT_VECTOR, globdat)
        globdat.model.takeAction(Action.GET_CONSTRAINTS, globdat)

        mbuild = globdat.get("mbuild")
        fext = globdat.get("fext")
        fint = globdat.get("fint")
        cons = globdat.get("cons")
        disp = globdat.get("solu")

        old_disp = deepcopy(disp)
        cons.updateSolution(disp)

        du = disp - old_disp
        K = mbuild.getDenseMatrix()
        r = fext - fint - np.array(K).dot(du)

        solver = Solver(self.type, cons)
        fdof = cons.getFdof()

        for iter in range(self.niter):

            # Update displacement vector
            solver.solve(K, du, r, mbuild.hbw)
            disp[fdof] += du[fdof]

            # Find interal force vector
            globdat.set("fint", np.zeros(ndof))

            if self.nrkey == "full":
                globdat.model.takeAction(Action.GET_MATRIX_0, globdat)
            elif self.nrkey == "mod" or self.nrkey == "LE":
                globdat.model.takeAction(Action.GET_INT_VECTOR, globdat)
            else:
                raise ValueError("{} not implemented !".format(self.nrkey))

            # Find out-of-balance force vector
            r = fext - globdat.get("fint")
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

    def __init__(self, name="sample"):
        self.name = name

    def init(self, conf, props, globdat):
        myProps = props.getProps(self.name)
        myConf = conf.makeProps(self.name)

        self.path = myProps.get("file")
        self.dofs = myProps.get("dofs")
        
        myConf.set("file", self.path)
        myConf.set("dofs", self.dofs)
        open(self.path, "w+").close()

        self.run(globdat)
        return Status.OK

    def run(self, globdat):
        i = globdat.i
        fint = globdat.get("fint")
        disp = globdat.get("solu")

        f = fint[self.dofs]
        u = disp[self.dofs]
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
    conf.writeFile("Examples/conf.pro")

    globdat.model.takeAction("PLOT_MESH", globdat)
    # mesh = globdat.get("mesh")
    # mesh.printDofSpace(15)
    # Run modules until Status.EXIT is issued
    while status != Status.EXIT:
        print("==== Running modules =====================")
        status = module.run(globdat)

    # Shutdown remaining modules
    print("==== Shutting down modules ===============")
    module.shutdown(globdat)
