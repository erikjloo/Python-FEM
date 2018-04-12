import scipy as np
from mesh import Mesh
from dofspace import DofSpace
from constraints import Constraints
from properties import Properties, ElementType
from algebra import MatrixBuilder
from shapes import Line2, Tri3, Quad4


path = '/home/erik/Documents/Python/FEM/square.msh'
mesh = Mesh()
mesh.readMesh(path)
nnod = mesh.nnodes
nele = mesh.nele

print(("Mesh read with {} nodes and {} elements").format(nnod,nele))

mesh.plotMesh()
dofs = DofSpace(nnod,2)
dofs.addTypes(["u","v"])
dofs.addDofs(range(nnod),["u","v"])

doftypes = ["u", "v"]
ndof = dofs.dofCount()

u = np.empty(ndof)

shape = Tri3()

props = Properties(nele,mesh.props)
props.addMaterial("Matrix", E = 29000, v = 0.3)

mbuild = MatrixBuilder(ndof)

# Iterate over elements assigned to model
for iele in range(nele):
    inods = mesh.getNodes(iele)
    coords = mesh.getCoords(inods)
    idofs = dofs.getDofIndices(inods,doftypes)
    # pylint: disable = unbalanced-tuple-unpacking
    [etype,mat] = props.getProperties(iele,["etype","mat"])

    # get nodal displacements

    # multiply w[ip] times thickness

    # element stiffness matrix
    ndof = len(idofs)
    kele = np.zeros((ndof,ndof))

    E, v = mat.E, mat.v
    la = v*E/((1+v)*(1-2*v))
    mu = E/(2*(1+v))
    D = np.array([[la+2*mu, la, 0], [la, la+2*mu, 0], [0, 0, mu]])

    # shape = newShape(etype, shape)
    [B, w] = shape.getBmatrix(coords)
    for ip in range(shape.nIP):


        # strain, eledisp = shape.getstrain( ip, ie)
        
        # Get tangent stiffness matrix:
        # material.update(stress, D, strain, ip)

        # Compute element stiffness matrix (kele):
        kele += (B[ip].transpose() @ D @ B[ip])*w[ip]
     
        # Compute the element force vector (fint):
        # fint += w[ip] * B[ip].transpose @ stress

    # Add k to the global stiffness matrix (Ksys):
    mbuild.addBlock(idofs, idofs, kele)

    # Add fint to the global force vector (Fint):
    # Fint[idofs] += fint

# output:
# u, strain, sigma, f_int, K

def newShape(etype):
    if etype == ElementType.line2:
        return Line2()
    elif etype == ElementType.tri3:
        return Tri3()
    elif etype == ElementType.quad4:
        return Quad4()
