# Import Standard Libraries
import scipy as np
from abc import ABCMeta, abstractmethod

# Import Local Libraries
from algebra import det, inv
from intSchemes import gauss_legendre, newton_cotes


#===========================================================================
#   ShapeFactory
#===========================================================================


def ShapeFactory(props):

    props = props.getProps("shape")
    type = props.get("type")
    scheme = props.get("scheme")

    if type == "Line2":
        print("Creating Line2 with", scheme, "quadrature")
        return Line2(scheme)

    elif type == "Tri3":
        print("Creating Tri3 with", scheme, "quadrature")
        return Tri3(scheme)

    elif type == "Quad4":
        print("Creating Quad4 with", scheme, "quadrature")
        return Tri3(scheme)


#===========================================================================
#   Shape
#===========================================================================


class Shape(metaclass=ABCMeta):
    """ Abstract Shape Class

    Virtual Static Members:
        nnod = number of nodes
        nIP = number of integration points
        ndim = number of dimensions of local coordinates

    Virtual Instance Members:
        N = list of array of shape functions at each IP
        N_xi = list of array of local shape gradients at each IP
        gp = array of integration points in local coordinates
        w = array of integration weights for each gp

    Virtual Public Methods:
        N(xi) = evalShapeFunctions(xi)
        N_xi(xi) = evalLocalGradients(xi)

    Public Methods:
        N[IP] = getShapeFunctions(IP=None)
        N_xi[IP] = getLocalGradients(IP=None)
        N_x[IP] = getGlobalGradients(coords, IP=None)
        [J,j] = getJacobian(coords, IP)
        x = getGlobalPoint(coords, xi)
        x[IP] = getGlobalIntegrationPoints(coords, IP=None)
        
    Private Methods:
        __setIntegrationScheme(scheme=None)
    """

    ndim = nIP = 0

    # Public:

    #-----------------------------------------------------------------------
    #   evalShapeFunctions
    #-----------------------------------------------------------------------

    @abstractmethod
    def evalShapeFunctions(self, xi):
        pass

    #-----------------------------------------------------------------------
    #   evalLocalGradients
    #-----------------------------------------------------------------------

    @abstractmethod
    def evalLocalGradients(self, xi):
        pass

    #-----------------------------------------------------------------------
    #   __init__
    #-----------------------------------------------------------------------

    def __init__(self, scheme=None):
        """ Input:  scheme = name of integration scheme
            Local:  N = list of [1 x nIP] arrays of shape functions at each IP
                    N_xi = list of [ndim x nIP] arrays of shape gradients at each IP """

        self.__setIntegrationScheme(scheme)
        # Shape functions and shape gradients
        self.N = []
        self.N_xi = []

        for ip in range(self.nIP):
            # Evaluate shape functions at each IP
            N = self.evalShapeFunctions(self.gp[ip])
            self.N.append(N)

            # Evaluate local shape gradients at each IP
            N_xi = self.evalLocalGradients(self.gp[ip])
            self.N_xi.append(N_xi)

    #-----------------------------------------------------------------------
    #   getShapeFunctions
    #-----------------------------------------------------------------------

    def getShapeFunctions(self, IP=None):
        """ Input:  IP = integration point number
            Output: N = array of shape functions at given IP """

        if IP is None:
            return self.N
        else:
            return self.N[IP]

    #-----------------------------------------------------------------------
    #   getLocalGradients
    #-----------------------------------------------------------------------

    def getLocalGradients(self, IP=None):
        """ Input:  IP = integration point number
            Output: N_xi = array of local shape gradients at given IP
                    w = weight at given IP """

        if IP is None:
            return [self.N_xi, self.w]
        else:
            return [self.N_xi[IP], self.w[IP]]

    #-----------------------------------------------------------------------
    #   getGlobalGradients
    #-----------------------------------------------------------------------

    def getGlobalGradients(self, coords, IP=None):
        """ Input:  coords = global coordinates of shape nodes
            Output: N_x = array of global shape gradients at given IP
                    w = w*j = weight at given IP """

        if coords.ndim > 1 and np.size(coords, 1) > self.ndim:
            raise ValueError("Element dimensions exceeded!")

        N_x = []
        w = []
        for ip in range(self.nIP):
            [J, j] = self.getJacobian(coords, ip)  # J = N_xi * coords
            DelN = np.dot(inv(J), self.N_xi[ip])  # DelN = inv(J) * N_xi
            w.append(self.w[ip]*j)
            N_x.append(DelN)

        if IP is None:
            return [N_x, w]
        else:
            return [N_x[IP], w[IP]]

    #-----------------------------------------------------------------------
    #   __init__
    #-----------------------------------------------------------------------

    def getIntegrationWeights(self, coords, IP=None):
        """ Input:  coords = global coordinates of shape nodes
            Output: w = w*j = weight at given IP """

        w = []
        for ip in range(self.nIP):
            [_, j] = self.getJacobian(coords, ip)
            w.append(self.w[ip]*j)

        if IP is None:
            return w
        else:
            return w[IP]

    #-----------------------------------------------------------------------
    #   getJacobian
    #-----------------------------------------------------------------------

    def getJacobian(self, coords, IP):
        """ Input:  coords = global coordinates of shape nodes
                    IP = integration point
            Output: J = N_xi * coords at given IP """

        J = self.N_xi[IP].dot(coords)
        j = det(J)
        if j < 0:
            raise ValueError("Negative jacobian!")
        return [J, j]

    #-----------------------------------------------------------------------
    #   getGlobalPoint
    #-----------------------------------------------------------------------

    def getGlobalPoint(self, xi, coords):
        """ Input:  xi = point in local coordinates
                    coords = global coordinates of shape nodes
            Output: x = array of global coordinates at given point """
        N = self.evalShapeFunctions(xi)
        x = N.dot(coords)
        return x

    #-----------------------------------------------------------------------
    #   getGlobalIntegrationPoints
    #-----------------------------------------------------------------------

    def getGlobalIntegrationPoints(self, coords, IP=None):
        """
        Input:
            coords = global coords of shape nodes
        Output:
            x = list of arrays of global coordinates of each IP
        """
        x = []
        for ip in range(self.nIP):
            x.append(self.N[ip].dot(coords))
        if IP is None:
            return x
        else:
            return x[IP]

    # Private:

    #-----------------------------------------------------------------------
    #   __setIntegrationScheme
    #-----------------------------------------------------------------------

    def __setIntegrationScheme(self, scheme=None):
        """
        Input:
            scheme = name of integration scheme
        Output:
            gp = array of weights and integration point local coords
        """
        if (scheme is None):  # Default is gauss_legendre
            [self.w, self.gp] = gauss_legendre(self.ndim, self.nIP)
        elif (scheme == "Gauss" or scheme == "GaussLegendre"):
            [self.w, self.gp] = gauss_legendre(self.ndim, self.nIP)
        elif (scheme == "Newton" or scheme == "NewtonCoates"):
            [self.w, self.gp] = newton_cotes(self.ndim, self.nIP)
        else:
            raise NotImplementedError('Int. Scheme does not exist.')


#===========================================================================
#   Line2
#===========================================================================


class Line2(Shape):
    """
    Line 2 element described by the following shape functions
        N1 = 0.5*(xi-1)
        N2 = 0.5*(xi+1)

    Static Members:
        nnod = number of nodes
        nIP = number of integration points
        ndim = number of dimensions of local coordinates

    Instance Members:
        N = list of array of shape functions at each IP
        N_xi = list of array of local shape gradients at each IP
        gp = array of integration points in local coordinates
        w = array of integration weights for each gp

    Public Methods:
        Line2()
        N(xi)= evalShapeFunctions(xi)
        N_xi(xi) = evalLocalGradients(xi)
        xi = getLocalPoint(x, coords)
    """

    # Static:
    nnod = 2  # number of nodes
    nIP = 1  # number of integration points
    ndim = 1  # number of dimensions (rank)

    # Public:
    def evalShapeFunctions(self, xi):
        """ Input: xi = point in local coordinates
            Output: [1 x 2] array of shape functions at given point """
        N1 = 0.5*(1-xi)
        N2 = 0.5*(1+xi)
        return np.array([N1, N2])

    def evalLocalGradients(self, xi=None):
        """ Output: [1 x 2] array of shape gradients at given point """
        return np.array([-0.5, 0.5])

    def getLocalPoint(self, x, coords):
        """ Input:  x = point in global coordinates
                    coord = global coordinates of shape nodes
            Output: xi = point in local coordinates """
        if coords.ndim is 1:
            return 2*(x - coords[0])/(coords[1]-coords[0]) - 1

        elif np.size(coords, 0) is 2:
            dx = coords[0, 1]-coords[0, 0]
            if np.abs(dx) < 1e-12:
                dy = coords[1, 1]-coords[1, 0]
                return 2*(x[1] - coords[1, 0])/dy - 1
            else:
                return 2*(x[0] - coords[0, 0])/dx - 1

        elif np.size(coords, 0) is 3:
            raise NotImplementedError("Not implemented for 3D")
        else:
            raise ValueError("Element dimensions exceeded!")


#===========================================================================
#   Line3
#===========================================================================


class Line3(Shape):
    """
    Line 3 element described by the following shape functions
        N1 = 0.5*xi^2-0.5*xi
        N2 = -xi^2+1
        N3 = 0.5*xi^2+0.5*xi

    Static Members:
        nnod = number of nodes
        nIP = number of integration points
        ndim = number of dimensions of local coordinates

    Instance Members:
        N = list of array of shape functions at each IP
        N_xi = list of array of local shape gradients at each IP
        gp = array of integration points in local coordinates
        w = array of integration weights for each gp

    Public Methods:
        Line3()
        N(xi)= evalShapeFunctions(xi)
        N_xi(xi) = evalLocalGradients(xi)
    """

    # Static:
    nnod = 3  # number of nodes
    nIP = 2  # number of integration points
    ndim = 1  # number of dimensions (rank)

    # Public:
    def evalShapeFunctions(self, xi):
        """ Input:  xi = point in local coordinates
            Output: [1 x 3] array of shape functions at given point """
        N1 = 0.5*xi**2-0.5*xi
        N2 = -xi**2+1
        N3 = 0.5*xi**2+0.5*xi
        return np.array([N1, N2, N3])

    def evalLocalGradients(self, xi):
        """ Input:  xi = point in local coordinates
            Output: [1 x 3] array of shape gradients at given point """
        N11 = xi-0.5
        N21 = -2*xi
        N31 = xi+0.5
        return np.array([N11, N21, N31])


#===========================================================================
#   Tri3
#===========================================================================


class Tri3(Shape):
    """
    Tri 3 element described by the following shape functions
        N1 = 1 - eta - xi
        N2 = xi
        N3 = eta

        3 o
          |\
          | \
          |  \
        1 o---o 2 --> xi
    
    Static Members:
        nnod = number of nodes
        nIP = number of integration points
        ndim = number of dimensions of local coordinates

    Instance Members:
        N = list of array of shape functions at each IP
        N_xi = list of array of local shape gradients at each IP
        gp = array of integration points in local coordinates
        w = array of integration weights for each gp

    Public Methods:
        Tri3()
        N(point)= evalShapeFunctions(point)
        N_xi(point) = evalLocalGradients(point)
        B = getBmatrix(coords)
    """

    # Static:
    nnod = 3
    nIP = 1
    ndim = 2

    # Public:
    def evalShapeFunctions(self, point):
        """ Input:  point = point in local coordinates
            Output: [1 x 3] array of shape functions at given point """
        xi = point[0]
        eta = point[1]
        N1 = 1 - eta - xi
        N2 = xi
        N3 = eta
        return np.array([N1, N2, N3])

    def evalLocalGradients(self, point=None):
        """ Output: [2 x 3] array of local shape gradients at given point """
        return np.array([[-1, 1, 0],
                         [-1, 0, 1]])

    def getBmatrix(self, coords):
        """ Input:  coords = global coords of shape nodes
            Output: list of [3 x 6] B matrices """
        [N_x, w] = self.getGlobalGradients(coords)
        B = []
        for IP in range(self.nIP):
            N = N_x[IP]
            B.append(np.array([[N[0, 0], 0, N[0, 1], 0, N[0, 2], 0],
                               [0, N[1, 0], 0, N[1, 1], 0, N[1, 2]],
                               [N[1, 0], N[0, 0], N[1, 1], N[0, 1], N[1, 2], N[0, 2]]]))
        return [B, w]


#===========================================================================
#   Quad4
#===========================================================================


class Quad4(Shape):
    """
    Quad 4 element described by the following shape functions
        N1 = 0.25*(1-xi-eta+xi*eta)
        N2 = 0.25*(1+xi-eta-xi*eta)
        N3 = 0.25*(1+xi+eta+xi*eta)
        N4 = 0.25*(1-xi+eta-xi*eta)

    Static Members:
        nnod = number of nodes
        nIP = number of integration points
        ndim = number of dimensions of local coordinates

    Instance Members:
        N = list of array of shape functions at each IP
        N_xi = list of array of local shape gradients at each IP
        gp = array of integration points in local coordinates
        w = array of integration weights for each gp

    Public Methods:
        Quad4()
        N(point)= evalShapeFunctions(point)
        N_xi(point) = evalLocalGradients(point)
        B = getBmatrix(coords)
    """
    
    # Static:
    nnod = 4  # number of nodes
    nIP = 4  # number of integration points
    ndim = 2  # number of dimensions

    # Public:
    def evalShapeFunctions(self, point):
        """ Input:  point = point in local coordinates
            Output: [1 x 4] array of shape functions at given point """
        xi = point[0]
        eta = point[1]
        N1 = 0.25*(1-xi-eta+xi*eta)
        N2 = 0.25*(1+xi-eta-xi*eta)
        N3 = 0.25*(1+xi+eta+xi*eta)
        N4 = 0.25*(1-xi+eta-xi*eta)
        return np.array([N1, N2, N3, N4])

    def evalLocalGradients(self, point):
        """ Input:  point = point in local coordinates
            Output: [2 x 4] array of local shape gradients at given point """
        xi = point[0]
        eta = point[1]
        N11 = 0.25*(-1+eta)
        N12 = 0.25*(-1+xi)
        N21 = 0.25*(1-eta)
        N22 = 0.25*(-1-xi)
        N31 = 0.25*(1+eta)
        N32 = 0.25*(1+xi)
        N41 = 0.25*(-1-eta)
        N42 = 0.25*(1-xi)
        return np.array([[N11, N21, N31, N41],
                         [N12, N22, N32, N42]])

    def getBmatrix(self, coords):
        """ Input:  coords = global coords of shape nodes
            Output: list of [3 x 8] B matrices """
        [N_x, w] = self.getGlobalGradients(coords)
        B = []
        for IP in range(self.nIP):
            N = N_x[IP]
            B.append(np.array([[N[0, 0], 0, N[0, 1], 0, N[0, 2], 0, N[0, 3], 0],
                               [0, N[1, 0], 0, N[1, 1], 0, N[1, 2], 0, N[1, 3]],
                               [N[1, 0], N[0, 0], N[1, 1], N[0, 1], N[1, 2], N[0, 2], N[1, 3], N[0, 3]]]))
        return [B, w]


#===========================================================================
#   Quad4
#===========================================================================
