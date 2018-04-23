# Import Standard Libraries
import numpy as np
from numpy import ix_
from numpy.linalg import norm

# Import Local Libraries
from solvers import Solver
from algebra import MatrixBuilder


def newton_raphson(model, cons, nrkey, niter, f_ext, f_int, u, K):
    """ Newton Raphson

    Directory of variables:
      Input:
        model
        fdof = vector of free degrees of freedom
        sdof = vector of supported (prescribed) degrees of freedom
        f_ext = external force vector of current load step
        f_int = internal force vector from previous load step
        u = displacement vector
        K = initial stiffness matrix
        niter = max. number of iterations per step
        nrkey = key for type of Newton-Raphson: 'full','mod','LE'
      Output:
        u = displacement vector at the end of the load step
        f_int = internal force vector at the end of the load step
        K = stiffness matrix for next load step
      Internal Variables:
        Da = change in displacement vector
        da = Iterative increment in displacement vector
        r = out-of-balance force vector
    """

    ndof = cons.dofCount()
    fdof = cons.get_fdof()

    # Initialize Data: Da = 0, r = f_ext - f_int
    Da = np.zeros(ndof)
    da = np.zeros(ndof)

    r = f_ext - f_int


    for i in range(niter):

        # Solve K_sys*da = r

        solver = Solver("rtfreechol", cons)

        da = solver.solve(K, da, r)

        # Update displacement vector
        Da[fdof] = Da[fdof] + da[fdof]

        # Update geometry
        # model.updateGeometry(d+Da)

        # Find internal strains

        # Find internal stresses

        # Find interal force vector

        # Find out-of-balance force vector
        r = f_ext - f_int
        nrm = norm(r[fdof])

        # Check convergence
        if i == 0:
            nrm1 = norm(r[fdof])

        if nrm < tol*nrm1:
            break

        # Update stiffness matrix
        if nrkey == 'full':
            mbuild = MatrixBuilder(model.dofCount())
            K = model.get_Matrix_0(mbuild, f_int)[0]

    u[fdof] = u[fdof] + Da[fdof]
    return(u,f_int,K)


def multistep(model, nsteps=1, nrkey="full", niter=20):
    """ Performs multi-step analysis of a given model

    Directory of variables:
      Input:

        nsteps = number of steps
        nrkey = key for type of Newton-Raphson: 'full','mod','LE'
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
    ndof = model.dofCount()

    # Apply Dirichlet boundary conditions
    d = model.get_Dirichlet()

    # Apply Neumann boundary conditions
    f = model.get_Neumann()

    # Matrix builder
    mbuild = MatrixBuilder(ndof)

    # Internal force vector
    f_int = np.zeros(ndof)

    # Displacement vector
    u = np.zeros(ndof)

    # Initial stiffness matrix
    K,f_int = model.get_Matrix_0(mbuild, f_int)

    # Allocate matrices for results
    F = []
    U = []

    # Loop over n steps
    for step in range(nsteps):
        f_ext = f*step/(nsteps+1)
        u[sdof] = d[sdof]*step/(nsteps+1)
        u, f_int, K = newton_raphson(model, nrkey, niter, fdof, sdof, f_ext, f_int, u, K)
        print(u)
        print(f_int)
        U = U.append(u.transpose())
        F = F.append(f_int.transpose())
    return(U,F)
