# Import Standard Libraries
import numpy as np
from copy import deepcopy
from numpy.linalg import norm

# Import Local Libraries
from solvers import Solver
from constraints import Constraints
from algebra import MatrixBuilder


def newton_raphson(disp, fext, cons, mbuild, fint, model, mesh, nrkey, niter, tol=1e-3):
    """ Newton Raphson

    Input:
        disp = displacement vector
        fext = external force vector of current load step
        mbuild = 
        fint = internal force vector from previous load step
        model =
        niter = max. number of iterations per step
        nrkey = key for type of Newton-Raphson: 'full','mod','LE'
    Output:
        u = displacement vector at the end of the load step
        fint = internal force vector at the end of the load step
        K = stiffness matrix for next load step
    Internal Variables:
        Da = change in displacement vector
        da = Iterative increment in displacement vector
        r = out-of-balance force vector
    """

    ndof = cons.dofCount()
    fdof = cons.get_fdof()

    # Initialize Data: Da = 0, r = fext - fint
    Da = np.zeros(ndof)
    da = np.zeros(ndof)

    r = fext - fint


    for i in range(niter):

        solver = Solver("rtfreechol", cons)
        K = mbuild.getDenseMatrix()
        da = solver.solve(K, da, r)

        # Update displacement vector
        Da[fdof] = Da[fdof] + da[fdof]
        # mesh.updateGeometry(d+Da)

        # Find interal force vector
        if nrkey == "full":
            model.get_Matrix_0(mbuild, fint, disp, mesh)
        elif nrkey == "mod" or nrkey == "LE":
            model.get_Int_Vector(fint, disp, mesh)

        # Find out-of-balance force vector
        r = fext - fint
        nrm = norm(r[fdof])

        # Check convergence
        if i == 0: nrm1 = deepcopy(nrm)
        if nrm < tol*nrm1: break

    disp[fdof] = disp[fdof] + Da[fdof]
    return [disp, fint]


def multistep(model, mesh, nsteps=1, nrkey="none", niter=20):
    """ Performs multi-step analysis of a given model

    Directory of variables:
      Input:
        nsteps = number of steps
        nrkey = key for type of Newton-Raphson: "none", "full", "mod", "LE"
        niter = max. number of iterations per step
      Output:
        U = array of displacement results for every load step
        F = array of force results for every load step
      Internal Variables:
        u = vector of displacements & rotations
        f = vector of forces & moments
        fdof = vector of free degrees of freedom
        sdof = vector of supported (prescribed) degrees of freedom
    """

    # Global Data
    ndof = mesh.dofCount()
    cons = Constraints(ndof)
    mbuild = MatrixBuilder(ndof)
    fint = np.zeros(ndof)
    fext = np.zeros(ndof)
    disp = np.zeros(ndof)

    # Constraints
    idofs = mesh.getDofIndices([1, 4], ["u", "v", "w"])
    cons.addConstraints(idofs)
    idofs = mesh.getDofIndices([2, 5], ["v", "w"])
    cons.addConstraints(idofs)
    idofs = mesh.getDofIndices([3, 6], ["u", "w"])
    cons.addConstraints(idofs)

    # Loads
    idofs = mesh.getDofIndices([3, 6], ["v"])
    fext[idofs] = 5

    # Initial stiffness matrix
    model.get_Matrix_0(mbuild, fint, disp, mesh)
    K = mbuild.getBlock(range(10),range(10))
    print(K)

    # Allocate matrices for results
    F = []
    U = []

    # Loop over n steps
    for step in range(nsteps):
        # fext = Load Model
        # disp = Disp Model
        newton_raphson(disp, fext, cons, mbuild, fint, model, mesh, nrkey, niter)
        U = U.append(disp)
        F = F.append(fint)
    return [U, F]
