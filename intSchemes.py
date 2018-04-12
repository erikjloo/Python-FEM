# Import Standard Libraries
import numpy as np


#===========================================================================
#   GaussLegendre
#===========================================================================


def gauss_legendre(dim,nIP):
    """
    Input:
        dim = number of dimensions of parent shape
        nIP = total number of integration points
    Output:
        w = weights per gauss point
        gp = local coordinates of gauss point
    """
    if (dim == 1):
        if (nIP == 1):	 # 1-point Gauss integration in 1D
            # Weights
            w = np.array([2])
            # Locations
            gp = np.array([0])
        elif (nIP == 2): # 2-point Gauss integration in 1D
            # Weights
            w = np.array([1, 1])
            # Locations
            a = 1/np.sqrt(3)
            gp = np.array([-a, a])
        elif (nIP == 3): # 3-point Gauss integration in 1D
            # Weights
            w = np.array([5/9,8/9,5/9])
            # Locations
            a = np.sqrt(3/5)
            gp = np.array([-a, 0, a])
        elif (nIP == 4): # 4-point Gauss integration in 1D
            w = np.zeros(4);             gp = np.zeros(4)
            w[0] = 0.6521451548625461;   gp[0] = -0.3399810435848563
            w[1] = 0.6521451548625461;   gp[1] = 0.3399810435848563
            w[2] = 0.3478548451374538;   gp[2] = -0.8611363115940526
            w[3] = 0.3478548451374538;   gp[3] = 0.8611363115940526
        else:
            raise ValueError('Integration scheme not programmed')
    elif (dim == 2):
        if (nIP == 1):	# 1-point Gauss integration in 2D for triangle
            # Weights
            w = np.array([0.5])
            # Locations
            gp = np.array([[1/3,1/3]])
        elif (nIP == 3): # 3-point Gauss integration in 2D for triangle
            # Weights
            w = np.array([1/6, 1/6, 1/6])
            # Locations
            gp = np.array([[0,0.5],
                           [0.5,0],
                           [0.5,0.5]])
        elif (nIP == 4): # 2x2-point Gauss integration in 2D for quad
            a = 1/np.sqrt(3)
            # Weights
            w = np.array([1, 1, 1, 1])
            # Locations
            gp = np.array([[-a,-a],
                           [a,-a],
                           [a,a],
                           [-a,a]])
        elif (nIP == 9): # 3x3-point Gauss integration in 2D for quad
            gp = np.zeros((3,9))
            a = np.sqrt(3/5)
            # Weights
            w = np.zeros(9)
            gp[0:3] = 25.0/81.0
            gp[4:7] = 40.0/81.0
            gp[8] = 64.0/81.0
            # Locations
            gp = np.array([[-a,-a],
                           [a,-a],
                           [a,a],
                           [-a,a],
                           [0,-a],
                           [a,0],
                           [0,a],
                           [-a,0],
                           [0,0]])
        else:
            raise NotImplementedError('Integration scheme not programmed')
    elif (dim == 3):
        if (nIP == 1): # 1-point Gauss integration in 3D for tetrahedron
            w = np.array([1/6])
            gp = np.array([[1/4, 1/4, 1/4]])
        else:
            raise NotImplementedError('Integration scheme not programmed')

    return w, gp


#===========================================================================
#   NewtonCotes
#===========================================================================


def newton_cotes(dim,nIP):
    """
    Input:
        dim = number of dimensions of parent shape
        nIP = total number of integration points
    Output:
        w = weights per gauss point
        gp = local coordinates of gauss point
    """
    if (dim == 1):
        if (nIP == 1):	 # 1-point Gauss integration in 1D
            # Weights
            w = np.array([2])
            # Locations
            gp = np.array([0])
        elif (nIP == 2): # 2-point Gauss integration in 1D
            # Weights
            w = np.array([1, 1])
            # Locations
            gp = np.array([-1, 1])
        elif (nIP == 3): # 3-point Gauss integration in 1D
            # Weights
            w = np.array([1/3,4/3,1/3])
            # Locations
            gp = np.array([-1, 0, 1])
        elif (nIP == 4): # 4-point Gauss integration in 1D
            w = np.array([1/4,3/4,3/4,1/4])
            gp = np.array([-1,-1/3,1/3,1])
        else:
            raise NotImplementedError('Integration scheme not programmed')
    elif (dim == 2):
        if (nIP == 1):	# 1-point Gauss integration in 2D for triangle
            # Weights
            w = np.array([0.5])
            # Locations
            gp = np.array([1/3,1/3])
        elif (nIP == 3):
            # Weights
            w = np.array([1/6, 1/6, 1/6])
            # Locations
            gp = np.array([[0, 0.5, 0.5],
                           [0.5, 0, 0.5]])
        elif (nIP == 4): # 2x2-point Gauss integration in 2D for quad
            a = 1/np.sqrt(3)
            # Weights
            w = np.array([1, 1, 1, 1])
            # Locations
            gp = np.array([[-a, a, a, -a],
                           [-a, -a, a, a]])
        elif (nIP == 9): # 3x3-point Gauss integration in 2D for quad
            gp = np.zeros((3,9))
            a = np.sqrt(3/5)
            # Weights
            w = np.zeros(9)
            gp[0:3] = 25.0/81.0
            gp[4:7] = 40.0/81.0
            gp[8] = 64.0/81.0
            # Locations
            gp = np.zeros((2,9))
            gp = np.array([[-a, a, a, -a, 0, a, 0, -a, 0],
                           [-a, -a, a, a, -a, 0, a, 0, 0]])
        else:
            raise NotImplementedError('Integration scheme not programmed')
    elif (dim == 3):
        if (nIP == 1): # 1-point Gauss integration in 3D for tetrahedron
            w = np.array([1/6])
            gp = np.array([[1/4, 1/4, 1/4]])
        else:
            raise NotImplementedError('Integration scheme not programmed')

    return w, gp
