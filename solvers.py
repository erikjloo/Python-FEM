# Import Standard Libraries
import numpy as np
import itertools
from numpy import ix_
from numpy.linalg import cond


#===========================================================================
#   Solver
#===========================================================================


class Solver(object):

    def __init__(self, method, cons):

        self.method = method
        self.ndof = cons.dofCount()
        self.fdof = cons.get_fdof()

    def solve(self, K, lhs, rhs, hbw=None):
        """ Solves Ax = b """
        hbw is self.ndof if hbw is None else hbw

        if self.method is "rtfreechol":
            lhs[self.fdof] = rtfreechol(K[ix_(self.fdof, self.fdof)], rhs[self.fdof], hbw)[0]
        elif self.method is "gauss_seidel":
            lhs[self.fdof] = gauss_seidel(K[ix_(self.fdof, self.fdof)], rhs[self.fdof])[0]

        return lhs, rhs

#===========================================================================
#   Root Free Cholesky
#===========================================================================


def rtfreechol(A,b,hbw):
    """ Solves Ax = b using root-free Cholesky

    Input:
        A = [n x n] symmetric matrix
        b = [1 x n] forcing vector
        hbw = half-band-width. Example:

            |-----|<--------hbw+1
            [x x x 0 0]
            [x x x x 0]
            [x x x x x]   --> hbw = 2
            [0 x x x x]
            [0 0 x x x]
            
    Local Variables:
        y = [1 x n] intermediate solution
        n = matrix size
    Output:
        x = [1 x n] solution vector
        L = [n x n] lower diagonal matrix
        D = [n x n] diagonal matrix
    """
    n = np.size(A,0)
    L = np.eye(n)       # create space for L
    D = np.eye(n)       # create space for D
    y = np.zeros(n)     # create space for intermediate solution y
    x = np.zeros(n)     # create space for solution vector x

    for j in range(n):

        D[j,j] = A[j,j]
        k_min = max([0,j-hbw])
        for k in range(k_min,j): # k = 0 : j
            D[j,j] -= D[k,k]*L[j,k]**2

        i_max = min([j+hbw+1,n])
        for i in range(j+1,i_max): # i = j+1 : n

            L[i,j] = A[i,j]
            k_min = max([0,i-hbw])
            for k in range(k_min,j): # k = 0 : j
                L[i,j] -= D[k,k]*L[i,k]*L[j,k]
            L[i,j] /= D[j,j]

        # Solve [L]{y} = {b}
        y[j] = b[j]
        i_min = max([0,j-hbw])
        for i in range(i_min,j):
            y[j] -= L[j,i]*y[i]
        y[j] /= L[j,j]

    
    for j in range(n-1,-1,-1):
        
        # Solve [D][L]'{x} = {y}
        x[j] = y[j]/D[j,j]
        i_max = min([j+hbw+1,n])
        for i in range(j+1,i_max):
            x[j] -= L[i,j]*x[i]
        x[j] /= L[j,j]

    if cond(A) > 10^18 or cond(A) < 10^-18:
        raise ValueError("Matrix A is ill-conditioned")

    return x,L,D


#===========================================================================
#   Gauss-Seidel
#===========================================================================


def gauss_seidel(A, b, beta=1, q=6, niter=200):
    """ Solves Ax = b using Gauss Seidel

    Input:
        A = [n x n] matrix with diagonal dominance
        b = [1 x n] forcing vector
        beta = relaxation or weighting factor
            1 < beta < 2 -- overrelaxation (speed up convergence rate)
            0 < beta < 1 -- underrelaxation (attempt to converge non-converging systems)
        q = desired digit accuracy
    Local:
        n = matrix size
        x_prev = prev. iteration of x
        x_new = estimate of x[i]
        E_max = maximum tolerance
        E_a = approximation error
    Output:
        x = [1 x n] solution vector
    """
    n = np.size(A,1)
    x = np.zeros(n)
    E_max = 0.5*10**(2-q)
    
    for _ in itertools.repeat(None, niter):
        
        x_prev = np.copy(x)
        
        for i in range(n):

            # Temporary x_new
            x_new = b[i]
            for k in range(0,i):
                x_new -= A[i,k]*x[k]
            for k in range(i+1,n):
                x_new -= A[i,k]*x[k]
            x_new /= A[i,i]

            # New x[i] as weighted avg. of x_new and x_prev[i]
            x[i] = beta*x_new + (1-beta)*x_prev[i]

        # Check convergence
        E_a = max(abs((x-x_prev)/x))*100
        if E_a < E_max: break
         
    return x


#===========================================================================
#   Example
#===========================================================================


if __name__ == "__main__":
    
    A = [[4.896,    0.0423,    0.6751,    0.0000,    0.0000],
        [0.0423,    9.7301,    0.1892,    0.9748,    0.0000],
        [0.6751,    0.1892,    6.2031,    0.6513,    0.2578],
        [0.0000,    0.9748,    0.6513,    2.3121,    0.4035],
        [0.0000,    0.0000,    0.2578,    0.4035,    1.5221]]

    A = np.array(A)*100

    b = [1,1,1,1,1]
    
    [x,L,D] = rtfreechol(A,b,2)
    print(x)
    x = np.linalg.solve(A,b)
    print(x)
    x = gauss_seidel(A,b)
    print(x)


