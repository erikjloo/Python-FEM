import scipy as np
from itemset import NodeSet,ElementSet
from dofspace import DofSpace
from shapes import Shape, Line2, Line3, Tri3, Quad4

example = "all"
if example == 1 or example == "all":
    #==Example 1=================================================#

    coords = np.array([[-1.2,-1],[1.4,-1],[1,1],[-1,1]])

    E,v = 1,0.25
    la = v*E/((1+v)*(1-2*v))
    mu = E/(2*(1+v))
    D = np.zeros((3,3))
    D[0,0] = D[1,1] = la+2*mu
    D[0,1] = D[1,0] = la
    D[2,2] = mu

    quad4 = Quad4()
    N = quad4.getShapeFunctions()
    IP = quad4.getGlobalIntegrationPoints(coords)
    [dN,w] = quad4.getGlobalGradients(coords)
    [B,w] = quad4.getBmatrix(coords)

    K = np.zeros((8,8))
    for ip in range(quad4.nIP):
        K += (B[ip].transpose() @ D @ B[ip])*w[ip]

    print("\n\n Quad4 \n IP =",IP)
    print(" Shape functions = \n",N)
    print(" Shape gradients = \n",dN)
    print(" weights = \n",w)
    print(" K  = \n",K)

if example == 2 or example == "all":
    #==Example 2=================================================#

    coords = np.array([[0,0],[1.2,0],[-0.1,1.3]])

    tri3 = Tri3("Gauss")
    N = tri3.getShapeFunctions()
    IP = tri3.getGlobalIntegrationPoints(coords)
    [dN,w] = tri3.getGlobalGradients(coords)
    [B,w] = tri3.getBmatrix(coords)
    K = np.zeros((6,6))
    for ip in range(tri3.nIP):
        K += B[ip].transpose() @ B[ip]*w[ip]

    print("\n\n Tri3 \n IP =",IP)
    print(" Shape functions = \n",N)
    print(" Shape gradients = \n",dN)
    print(" weights = \n",w)
    print(" K  = \n",K)

if example == 3 or example == "all":
    #==Example 3=================================================#

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
    #==Example 4=================================================#

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

    for e in range(nele):
        connect = elems2.getNodes(e)
        coords = nodes2.getCoords(connect)
        quad4 = Quad4("GaussLegendre")
        [B,w] = quad4.getGlobalGradients(coords)
        IP = quad4.getGlobalIntegrationPoints(coords)
        for ip in range(4):
            A += w[ip]*(IP[ip][0]**2+IP[ip][1]**2)

    print(A)
