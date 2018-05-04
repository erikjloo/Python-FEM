import scipy as np
from itemset import NodeSet,ElementSet
from dofspace import DofSpace
from shapes import Shape, Line2, Line3, Tri3, Quad4

example = 1
if example == 1 or example == "all":
    #==Example 1=================================================#
    """ Quad 4 Shape Stiffness Matrix """

    coords = np.array([[-1.2,-1],[1.4,-1],[1,1],[-1,1]])

    E,v = 1,0.25
    la = v*E/((1+v)*(1-2*v))
    mu = E/(2*(1+v))
    D = np.zeros((3,3))
    D[0,0] = D[1,1] = la+2*mu
    D[0,1] = D[1,0] = la
    D[2,2] = mu
    quad4 = Quad4()
    N = quad4.getNmatrix()
    IP = quad4.getGlobalIntegrationPoints(coords)
    [dN,w] = quad4.getGlobalGradients(coords)
    [B,w] = quad4.getBmatrix(coords)

    K = np.zeros((8,8))
    for ip in range(quad4.nIP):
        K += (B[ip].transpose() @ D @ B[ip])*w[ip]

    print("\n\n Quad4 \n")
    print(" D = \n",D)
    print(" Integration points = \n", IP)
    print(" Shape functions = \n",N)
    print(" Shape gradients = \n",dN)
    print(" weights = \n",w)
    print(" K  = \n",K)

if example == 2 or example == "all":
    #==Example 2=================================================#
    """ Tri 3 Shape Stiffness Matrix """

    coords = np.array([[0,0],[1.2,0],[-0.1,1.3]])

    tri3 = Tri3(scheme="Gauss")
    N = tri3.getShapeFunctions()
    IP = tri3.getGlobalIntegrationPoints(coords)
    [dN,w] = tri3.getGlobalGradients(coords)
    [B,w] = tri3.getBmatrix(coords)
    K = np.zeros((6,6))
    for ip in range(tri3.nIP):
        K += B[ip].transpose() @ B[ip]*w[ip]

    print("\n\n Tri3 \n")
    print(" Integration points = \n", IP)
    print(" Shape functions = \n",N)
    print(" Shape gradients = \n",dN)
    print(" weights = \n",w)
    print(" K  = \n",K)

if example == 3 or example == "all":
    #==Example 3=================================================#
    """ Line 3 Shape Stiffness Matrix """

    coords = np.array([2,5,8])

    line3 = Line3()
    N = line3.getShapeFunctions()
    IP = line3.getGlobalIntegrationPoints(coords)
    [dN,w] = line3.getGlobalGradients(coords)
    K = np.zeros((3,3))
    for ip in range(line3.nIP):
        K += w[ip]*np.outer( dN[ip], dN[ip] )

    print("\n\n Line3 \n Int. Point = \n",IP)
    print(" Shape functions = \n",N)
    print(" Shape gradients = \n",dN)
    print(" weights = \n",w)
    print(" K  = \n",K)

if example == 4 or example == "all":
    #==Example 5=================================================#
    """ Line 2 Boundary Element """

    coords = np.array([[0.2, 0.3, 0.5],
                       [0.8, 0.9, 0.5]])

    x = np.array([0.5, 0.7, 0.5])

    line2 = Line2(scheme="Gauss2")
    N = line2.getShapeFunctions()
    IP = line2.getGlobalIntegrationPoints(coords)
    w = line2.getIntegrationWeights(coords)
    xi = line2.getLocalPoint(x, coords)

    print("\n\n Line3 \n Int. Point = \n", IP)
    print(" Shape functions = \n", N)
    print(" weights = \n", w)
    print(" xi = \n", xi)

if example == 5 or example == "all":
    #==Example 4=================================================#
    """ Area of a Mesh """

    nodes2 = NodeSet()
    elems2 = ElementSet()

    nele, nv, nu = 4,2,5
    dx,dy = 1.2,1

    for ni in range(nv):
        for nj in range(nu):
            coord = [dx*nj,dy*ni]
            nodes2.addNode(coord)

    for e in range(nele):
        connect = [e,e+1,e+6,e+5]
        elems2.addElement(connect)

    A = 0

    quad4 = Quad4("GaussLegendre")
    for e in range(nele):
        connect = elems2.getNodes(e)
        coords = nodes2.getCoords(connect)
        [B,w] = quad4.getGlobalGradients(coords)
        IP = quad4.getGlobalIntegrationPoints(coords)
        for ip in range(4):
            A += w[ip]*(IP[ip][0]**2+IP[ip][1]**2)

    print("Area of the mesh is", A)
