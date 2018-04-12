import numpy as np
import solvers as ls

from numpy import ix_
from numpy.linalg import norm

def transformation_matrix(coord_ele,webdir_ele):
    """ Assembles transformation matrix gamma

    Directory of variables:
      Input:
        coord_ele = [[x_coord,y_coord],  <-- start
                     [x_coord,y_coord]]  <-- finish
        webdir_ele = [x_glo,y_glo]
      Ouput:
        gamm = 4 x 4 transfomation matrix for 2-noded axial element
    """
    gamma = np.zeros(4,4)
    dx_dy = coord_ele[1,:] - coord_ele[1,:]
    dx_dy = norm(dx_dy)
    rotation = np.array([[0,-1],[1,0]])
    webdir_ele = dot(rotation,dx_dy)
    gamma[0,0:2] = dx_dy/L
    gamma[1,0:2] = webdir_ele
    gamma[ix_([2,3],[2,3])] = gamma[ix_([0,1],[0,1])]
    return(gamma)

def newton_raphson(nnodes,coord,nele,ele,u_ext,f_ext,f_int,fdof,sdof,K,hbw,ndof,tol,niter,nrkey):
    """ Newton Raphson

    Directory of variables:
      Input:
        nnodes = number of nodes
        coord = initial values of x, y, z, thx, thy, & thz
        nele = number of elements
        ele = element end nodes / connectivity matrix
        u_ext = initial displacement matrix
        f_ext = external force vector of current load step
        f_int = internal force vector from previous load step
        fdof = vector of free degrees of freedom
        sdof = vector of supported (prescribed) degrees of freedom
        K = initial stiffness matrix
        hbw = half-band width
        ndof = number of degrees of freedom per node
        tol = convergence tolerance
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
    # Initialize Data: Da = 0, r = f_ext - f_int
    Da = np.zeros((nnodes*ndof,1))
    da = np.zeros((nnodes*ndof,1))
    r = f_ext - f_int

    for i in range(niter):

        # Solve K_sys*da = r
        da[fdof] = ls.rtfreechol(K[ix_(fdof,fdof)],r[fdof],hbw)[0]

        # Update displacement vector
        Da[fdof] = Da[fdof] + da[fdof]

        # Update geometry
        coord = coord + np.reshape(d+Da,[nnodes,ndof])

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
            K = assembly(nnodes,coord,nele,ele,ndof)[0]

    u[fdof] = u[fdof] + Da[fdof]
    return(u,f_int,K)


def multistep(nnodes,coord,nele,ele,disp,load,tol,niter,nrkey,nsteps):
    """ Performs multi-step analysis of a given frame

    Directory of variables:
      Input:
        nnodes = number of nodes
        coord = initial values of x, y, z, thx, thy, & thz
        nele = number of elements
        ele = element end nodes / connectivity matrix
        disp = array of prescribed degrees of freedom
        load = array of loads & moments
        tol = convergence tolerance
        niter = max. number of iterations per step
        nrkey = key for type of Newton-Raphson: 'full','mod','LE'
        nsteps = number of steps
      Output:
        U = array of displacement results for every load step
        F = array of force results for every load step
      Internal Variables:
        u = vector of displacements & rotations
        f = vector of forces & moments
        fdof = vector of free degrees of freedom
        sdof = vector of supported (prescribed) degrees of freedom
    """
    ndof = 2

    # Create vectors for solving f = K*d
    d = np.isnan(ndof)
    f = np.isnan(ndof)

    # Apply Dirichlet boundary conditions

    # Apply Neumann boundary conditions
    
    # Indices for free & supported dof
    fdof = np.argwhere(np.isnan(d))[:,0]
    sdof = np.argwhere(np.isfinite(d))[:,0]

    # Internal force vector for first load step
    f_int = np.zeros((ndof*nnodes,1))

    # Displacement vector for first load step
    u = np.zeros((ndof*nnodes,1))

    # Initial stiffness matrix
    K,hbw = assembly(nnodes,coord,nele,ele,ndof)

    # Allocate matrices for results
    U = []
    F = []

    # Loop over n steps
    for step in range(1,nsteps+1):
        u[sdof] = d[sdof]*step/nsteps
        f_ext = f*step/nsteps
        u,f_int = nr_full(nnodes,coord,nele,ele,u,f_ext,f_int,fdof,sdof,K,hbw,ndof,tol,niter)
        print(u)
        print(f_int)
        U = U.append(u.transpose())
        F = F.append(f_int.transpose())
    return(U,F)
