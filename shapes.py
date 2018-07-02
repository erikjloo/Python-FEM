# Import Standard Libraries
import re
import scipy as np
from warnings import warn
from abc import ABCMeta, abstractmethod

# Import Local Libraries
from algebra import determinant, inverse, norm, gram_schmidt
from intSchemes import gauss_legendre, newton_cotes

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
        N = list of arrays of shape functions at each IP
        N_xi = list of arrays of local shape gradients at each IP
        w = list of integration weights at each IP
        gp = list of local coordinates of each IP

    Virtual Public Methods:
        N(xi) = evalShapeFunctions(xi)
        N_xi(xi) = evalLocalGradients(xi)

    Public Methods:
        N = getShapeFunctions(IP=None)
        N = getNmatrix(IP=None)
        N(xi) = evalNmatrix(xi)
        [N_xi, w] = getLocalGradients(IP=None)
        [N_x, w] = getGlobalGradients(coords, IP=None)
        [B, w] = getBmatrix(coords, IP=None)
        [J, j] = getJacobian(coords, IP)
        x = getGlobalPoints(coords, IP=None)
        x(xi) = evalGlobalPoint(coords, xi)
        coords = getLocalCoords(coords)
        
    Private Methods:
        __setIntegrationScheme(scheme=None)
    """

    nIP = nnod = ndim = 0

    # Public:

    #-----------------------------------------------------------------------
    #   evalShapeFunctions
    #-----------------------------------------------------------------------

    @abstractmethod
    def evalShapeFunctions(self, xi):
        raise NotImplementedError()

    #-----------------------------------------------------------------------
    #   evalLocalGradients
    #-----------------------------------------------------------------------

    @abstractmethod
    def evalLocalGradients(self, xi):
        raise NotImplementedError()

    #-----------------------------------------------------------------------
    #   __init__
    #-----------------------------------------------------------------------

    def __init__(self, scheme="Gauss"):
        """ Input:  scheme = name of integration scheme
            Local:  N = list of arrays of shape functions at each IP
                    N_xi = list of arrays of shape gradients at each IP """

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
        """ Input:  IP = integration point
            Output: N = array of shape functions at given IP """

        if IP is None:  # ==================================================
            return self.N
        else:  # ===========================================================
            return self.N[IP]

    #-----------------------------------------------------------------------
    #   getLocalGradients
    #-----------------------------------------------------------------------

    def getLocalGradients(self, IP=None):
        """ Input:  IP = integration point
            Output: N_xi = array of shape gradients at given IP
                    w = weight at given IP """

        if IP is None:  # ==================================================
            return [self.N_xi, self.w]
        else:  # ===========================================================
            return [self.N_xi[IP], self.w[IP]]

    #-----------------------------------------------------------------------
    #   getNmatrix
    #-----------------------------------------------------------------------

    def getNmatrix(self, IP=None, ndim=None):
        """ Input:  IP = integration point
            Output: N = N matrix at given IP """

        n = self.getShapeFunctions(IP)

        if IP is None:  # ==================================================
            N = []
            if ndim is None:
                for IP in range(self.nIP):
                    H = np.eye(self.ndim)*n[IP][0]
                    for nod in range(1, self.nnod):
                        H = np.hstack((H, np.eye(self.ndim)*n[IP][nod]))
                    N.append(H)
            else:
                for IP in range(self.nIP):
                    H = np.eye(ndim)*n[IP][0]
                    for nod in range(1, self.nnod):
                        H = np.hstack((H, np.eye(ndim)*n[IP][nod]))
                    N.append(H)

        else:  # ===========================================================
            if ndim is None:
                N = np.eye(self.ndim)*n[0]
                for nod in range(1, self.nnod):
                    N = np.hstack((N, np.eye(self.ndim)*n[nod]))
            else:
                N = np.eye(ndim)*n[0]
                for nod in range(1, self.nnod):
                    N = np.hstack((N, np.eye(ndim)*n[nod]))

        return N

    #-----------------------------------------------------------------------
    #   evalNmatrix
    #-----------------------------------------------------------------------

    def evalNmatrix(self, xi, ndim=None):
        """ Input: xi = point in local coordinates
            Output: N = N matrix at given point """

        n = self.evalShapeFunctions(xi)

        if ndim is None:  # ================================================
            N = np.eye(self.ndim)*n[0]
            for nod in range(1, self.nnod):
                N = np.hstack((N, np.eye(self.ndim)*n[nod]))
        else:  # ===========================================================
            N = np.eye(ndim)*n[0]
            for nod in range(1, self.nnod):
                N = np.hstack((N, np.eye(ndim)*n[nod]))

        return N

    #-----------------------------------------------------------------------
    #   getGlobalGradients
    #-----------------------------------------------------------------------

    def getGlobalGradients(self, coords, IP=None):
        """ Input:  coords = global coordinates of shape nodes
                    IP = integration point
            Output: N_x = array of global shape gradients at given IP
                    w = w*j = weight at given IP """

        if coords.ndim > 1:
            if np.size(coords, axis=1) > self.ndim:
                raise ValueError("Element dimensions exceeded!")

        if IP is None:  # ==================================================
            N_x = []
            w = []
            for ip in range(self.nIP):
                [J, j] = self.getJacobian(coords, ip)  # J = N_xi * coords
                N_x.append( inverse(J).dot(self.N_xi[ip]))
                w.append(self.w[ip]*j)

        else:  # ===========================================================
            [J, j] = self.getJacobian(coords, IP)
            N_x = inverse(J).dot(self.N_xi[IP])
            w = self.w[IP]*j

        return [N_x, w]

    #-----------------------------------------------------------------------
    #   getBmatrix
    #-----------------------------------------------------------------------

    def getBmatrix(self, coords, IP=None):
        """ Input:  coords = global coordinates of shape nodes
                    IP = integration point
            Output: B = B (or dN) matrix at given IP
                    w = w*j = weight at given IP """
        [N_x, w] = self.getGlobalGradients(coords, IP)

        if self.ndim == 1:
            return [N_x, w]

        elif self.ndim == 2:

            if IP is None:  # ==============================================
                B = []
                for ip in range(self.nIP):
                    # Assemble B matrix at given IP
                    b = np.zeros((3, 2*self.nnod))
                    for nod in range(self.nnod):
                        b[0, 2*nod] = b[2, 2*nod+1] = N_x[ip][0, nod]
                        b[1, 2*nod+1] = b[2, 2*nod] = N_x[ip][1, nod]
                    B.append(b)

            else:  # =======================================================
                B = np.zeros((3, 2*self.nnod))
                for nod in range(self.nnod):
                    B[0, 2*nod] = B[2, 2*nod+1] = N_x[0, nod]
                    B[1, 2*nod+1] = B[2, 2*nod] = N_x[1, nod]

            return [B, w]

        elif self.ndim == 3:

            if IP is None:  # ==============================================
                B = []
                for ip in range(self.nIP):
                    # Assemble B matrix at given IP
                    b = np.zeros((6, 3*self.nnod))
                    for nod in range(self.nnod):
                        b[0, 3*nod] = b[3, 3*nod+1] = b[5,
                                                        3*nod+2] = N_x[ip][0, nod]
                        b[1, 3*nod+1] = b[3, 3*nod] = b[4,
                                                        3*nod+2] = N_x[ip][1, nod]
                        b[2, 3*nod+2] = b[4, 3*nod+1] = b[5,
                                                          3*nod] = N_x[ip][2, nod]
                    B.append(b)

            else:  # =======================================================
                B = np.zeros((6, 3*self.nnod))
                for nod in range(self.nnod):
                    B[0, 3*nod] = B[3, 3*nod+1] = B[5, 3*nod+2] = N_x[0, nod]
                    B[1, 3*nod+1] = B[3, 3*nod] = B[4, 3*nod+2] = N_x[1, nod]
                    B[2, 3*nod+2] = B[4, 3*nod+1] = B[5, 3*nod] = N_x[2, nod]

            return [B, w]

    #-----------------------------------------------------------------------
    #   getStrain
    #-----------------------------------------------------------------------

    def getStrain(self, coords, disp, IP=None):
        """ Input:  coords = global coordinates of shape nodes
                    disp = displacement vector of shape nodes
            Output: strain = strain vector at given IP
                    B = B (or dN) matrix at given IP
                    w = w*j = weight at given IP """

        [B, w] = self.getBmatrix(coords, IP)

        if IP is None:  # ==================================================
            strain = []
            for ip in range(self.nIP):
                strain.append(B[ip].dot(disp))
            return [strain, B, w]

        else:  # ===========================================================
            strain = B.dot(disp)
            return [strain, B, w]

    #-----------------------------------------------------------------------
    #   getIntegrationWeights
    #----------------------------------------------------------s-------------

    def getIntegrationWeights(self, coords, IP=None):
        """ Input:  coords = global coordinates of shape nodes
                    IP = integration point
            Output: w = w*j = weight at given IP """

        if coords.ndim > 1:
            if np.size(coords, axis=1) > self.ndim:
                coords = self.getLocalCoords(coords)

        if IP is None:  # ==================================================
            w = []
            for ip in range(self.nIP):
                [_, j] = self.getJacobian(coords, ip)
                w.append(self.w[ip]*j)

        else:  # ===========================================================
            [_, j] = self.getJacobian(coords, IP)
            w = self.w[IP]*j

        return w

    #-----------------------------------------------------------------------
    #   getJacobian
    #-----------------------------------------------------------------------

    def getJacobian(self, coords, IP):
        """ Input:  coords = global coordinates of shape nodes
                    IP = integration point
            Output: J = N_xi * coords at given IP 
                    j = determinant(J) at given IP """

        J = self.N_xi[IP].dot(coords)
        j = determinant(J)
        if j < 0:
            raise ValueError("Negative jacobian!")
        return [J, j]

    #-----------------------------------------------------------------------
    #   getGlobalIntegrationPoints
    #-----------------------------------------------------------------------

    def getGlobalPoints(self, coords, IP=None):
        """ Input:  coords = global coordinates of shape nodes
                    IP = integration point
            Output: x = global coordinates of given IP """

        if IP is None:  # ==================================================
            x = []
            for ip in range(self.nIP):
                x.append(self.N[ip].dot(coords))

        else:  # ===========================================================
            x = self.N[IP].dot(coords)

        return x

    #-----------------------------------------------------------------------
    #   getGlobalPoint
    #-----------------------------------------------------------------------

    def evalGlobalPoint(self, xi, coords):
        """ Input:  coords = global coordinates of shape nodes
                    xi = point in local coordinates
            Output: x = global coordinates of given point """

        N = self.evalShapeFunctions(xi)
        return N.dot(coords)

    # Private:

    #-----------------------------------------------------------------------
    #   __setIntegrationScheme
    #-----------------------------------------------------------------------

    def __setIntegrationScheme(self, scheme):
        """ Input:  scheme = name of integration scheme
            Local:  w = list of integration weights at each IP
                    gp = list of local coordinates of each IP """

        try:  # User may input scheme = "Gauss4", for example.
            r = re.compile("([a-zA-Z]+)([0-9]+)")
            m = r.match(scheme)
            scheme = m.group(1)
            self.nIP = int(m.group(2))
        except:
            AttributeError()

        if (scheme == "Gauss" or scheme == "GaussLegendre"):
            [self.w, self.gp] = gauss_legendre(self.ndim, self.nIP)
        elif (scheme == "Newton" or scheme == "NewtonCoates"):
            [self.w, self.gp] = newton_cotes(self.ndim, self.nIP)
        else:
            raise NotImplementedError('Int. Scheme does not exist.')

    #-----------------------------------------------------------------------
    #   getLocalCoords
    #-----------------------------------------------------------------------

    def getLocalCoords(self, coords):
        """ Input:  coords = matrix of global coordinates
            Output: coords = matrix of coordinates in local coordinates """

        if self.ndim == 1 and np.size(coords, axis=1) == 2:

            # Transformation matrix
            dx_dy = coords[1, :] - coords[0, :]
            i_bar = dx_dy/norm(dx_dy)
            j_bar = np.dot(np.array([[0, -1], [1, 0]]), i_bar)
            Gamma = [i_bar, j_bar]

            # Transform into local coordinates
            coords = Gamma @ coords.transpose()
            coords = coords.transpose()

            # Remove y-coordinates
            return coords[:, 0]

        elif self.ndim == 1 and np.size(coords, axis=1) == 3:

            # Transformation matrix
            i_bar = coords[1, :] - coords[0, :]  # [dx, dy, dz]
            j_bar = [0, 1, 0]  # Assume j-bar points upwards
            k_bar = np.cross(i_bar, j_bar)
            Gamma = gram_schmidt(i_bar, j_bar, k_bar)

            # Transform into local coordinates
            coords = Gamma @ coords.transpose()
            coords = coords.transpose()

            # Remove y and z-coordinates
            return coords[:, 0]

        elif self.ndim == 2 and np.size(coords, axis=1) == 3:

            warn("Unverified result, please verify!!")

            # Transformation matrix
            i_bar = coords[2, :] - coords[1, :]  # [dx, dy, dz along one side]
            # [dx, dy, dz along another side]
            j_bar = coords[3, :] - coords[2, :]
            k_bar = np.cross(i_bar, j_bar)
            i_bar, j_bar, k_bar = gram_schmidt(i_bar, j_bar, k_bar)
            Gamma = [i_bar, j_bar, k_bar]

            # Transform into local coordinates
            coords = Gamma @ coords.transpose()
            coords = coords.transpose()

            # Remove z-coordinates
            return coords[:, 0:2]

        elif self.ndim == 3 and np.size(coords, axis=1) > 3:
            raise ValueError("Element dimensions exceeded !")

    @staticmethod
    def shapeFactory(conf, props):
        myProps = props.getProps("shape")
        myConf = conf.makeProps("shape")

        type = myProps.get("type")
        scheme = myProps.get("scheme", "Gauss")

        myConf.set("type", type)
        myConf.set("scheme", scheme)

        message = "Creating {} with {} quadrature"
        print(message.format(type, scheme))

        if type == "Line2":
            return Line2(scheme)
        elif type == "Tri3":
            return Tri3(scheme)
        elif type == "Quad4":
            return Quad4(scheme)
        elif type == "Tetra4":
            return Tetra4(scheme)
        else:
            raise KeyError("{} shape not implemented".format(type))

#===========================================================================
#   Line2
#===========================================================================


class Line2(Shape):
    """
    Line 2 element described by the following shape functions
        N1 = 0.5*(xi-1)
        N2 = 0.5*(xi+1)

    1-------------2--> xi
    
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
    nIP = 1  # number of integration points
    nnod = 2  # number of nodes
    ndim = 1  # number of dimensions (localrank)

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
        if coords.ndim == 1:
            dx = coords[1] - coords[0]
            return 2*(x - coords[0])/dx - 1

        elif np.size(coords, axis=1) == 2 or np.size(coords, axis=1) == 3:

            # Find point on line closest to x
            x0 = coords[0]
            dx = coords[-1] - coords[0]
            t = -(x0.dot(dx) - x.dot(dx))/dx.dot(dx)

            # Map x onto line:
            x = x0 + t*dx

            # Transform to local coordinates
            coords = self.getLocalCoords(np.vstack((coords, x)))
            return self.getLocalPoint(coords[-1], coords[:-1])

        else:
            raise ValueError("Element dimensions exceeded !")


#===========================================================================
#   Line3
#===========================================================================


class Line3(Shape):
    """
    Line 3 element described by the following shape functions
        N1 = 0.5*xi^2-0.5*xi
        N2 = -xi^2+1
        N3 = 0.5*xi^2+0.5*xi

    1-------2------3--> xi

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
    nIP = 2  # number of integration points
    nnod = 3  # number of nodes
    ndim = 1  # number of dimensions (localrank)

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
        
          ^ eta
          |
          3
          |`| 
          |  `| 
          |    `| 
          |      `| 
          1--------2 --> xi
    
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
    """

    # Static:
    nIP = 1  # number of integration points
    nnod = 3  # number of nodes
    ndim = 2  # number of dimensions (localrank)

    # Public:
    def evalShapeFunctions(self, point):
        """ Input:  point = point in local coordinates
            Output: [1 x 3] array of shape functions at given point """
        xi, eta = point
        return np.array([1 - eta - xi, xi, eta])

    def evalLocalGradients(self, point=None):
        """ Output: [2 x 3] array of local shape gradients at given point """
        return np.array([[-1, 1, 0],
                         [-1, 0, 1]])


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

              eta
              ^
              |
        4-----------3 
        |     |     |  
        |     |     |  
        |     +---- | --> xi
        |           |  
        |           |  
        1-----------2  

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
    """

    # Static:
    nIP = 4  # number of integration points
    nnod = 4  # number of nodes
    ndim = 2  # number of dimensions (localrank)

    # Public:
    def evalShapeFunctions(self, point):
        """ Input:  point = point in local coordinates
            Output: [1 x 4] array of shape functions at given point """
        xi, eta = point
        N1 = 0.25*(1-xi-eta+xi*eta)
        N2 = 0.25*(1+xi-eta-xi*eta)
        N3 = 0.25*(1+xi+eta+xi*eta)
        N4 = 0.25*(1-xi+eta-xi*eta)
        return np.array([N1, N2, N3, N4])

    def evalLocalGradients(self, point):
        """ Input:  point = point in local coordinates
            Output: [2 x 4] array of local shape gradients at given point """
        xi, eta = point
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


#===========================================================================
#   Tetra4
#===========================================================================


class Tetra4(Shape):
    """
    Quad 4 element described by the following shape functions
        N1 = 1-xi-eta-zeta
        N2 = xi
        N3 = eta
        N4 = zeta


                   xi
                 .
               ,/
              /
           3                 
         ,/|`| 
       ,/  |  `| 
     ,/    '.   `|
   ,/       |     `|
 ,/         |       `|
1-----------'.--------2 --> eta
 `|.         |      ,/ 
    `|.      |    ,/
       `|.   '. ,/
          `|. |/
             `4 
                `|.
                   ` zeta
                   
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
        Tetra4()
        N(point)= evalShapeFunctions(point)
        N_xi(point) = evalLocalGradients(point)
    """

    # Static:
    nIP = 1  # number of integration points
    nnod = 4  # number of nodes
    ndim = 3  # number of dimensions (localrank)

    # Public:
    def evalShapeFunctions(self, point):
        """ Input:  point = point in local coordinates
            Output: [1 x 4] array of shape functions at given point """
        xi, eta, zeta = point
        return np.array([1-xi-eta-zeta, xi, eta, zeta])

    def evalLocalGradients(self, point):
        """ Input:  point = point in local coordinates
            Output: [3 x 4] array of local shape gradients at given point """
        return np.array([[-1, 1, 0, 0],
                         [-1, 0, 1, 0],
                         [-1, 0, 0, 1]])
