# Import Standard Libraries
import matplotlib.pyplot as plt
import scipy as np
import time

# Import Local Libraries
from properties import Properties
from modules import InputModule, InitModule
from solvers import Solver

np.set_printoptions(precision=4)

start = time.time()

dim = 3
# Initialization
if dim == 2:
    file = "Examples/2D_semicircle.pro"
elif dim == 3:
    file = "Examples/3D_semicircle.pro"

# Props
props = Properties()
props.parseFile(file)
props.print()

# Mesh
module = InputModule("input")
mesh = module.init(props)

# Model & Global Data
module = InitModule()
[model, cons, mbuild, fint, fext, disp] = module.init(props, mesh)

# Constraints
if dim == 2:
    idofs = mesh.getDofIndices(0, ["u", "v"])
    cons.addConstraints(idofs)
    idof = mesh.getDofIndex(2, "v")
    cons.addConstraint(idof)
    idof = mesh.getDofIndex(1, "u")
    cons.addConstraint(idof)

    # Loads
    idof = mesh.getDofIndex(1, "v")
    fext[idof] = 10

elif dim == 3:
    idofs = mesh.getDofIndices([0, 3], ["u", "v", "w"])
    cons.addConstraints(idofs)
    idofs = mesh.getDofIndices([1, 4], ["v", "w"])
    cons.addConstraints(idofs)
    idofs = mesh.getDofIndices([2, 5], ["u", "w"])
    cons.addConstraints(idofs)

    # Loads
    idofs = mesh.getDofIndices([2, 5], "v")
    fext[idofs] = 5

# Initial stiffness matrix
hbw = model.get_Matrix_0(mbuild, fint, disp, mesh)
print("The half-band-width is",hbw)

# Solve
solver = Solver("numpy", cons)
K = mbuild.getDenseMatrix()
solver.solve(K, disp, fext, hbw)

# K = mbuild.getMatrix()
# plt.spy(K)
# plt.show()

stop = time.time()
print("Elapsed time is ",stop-start)

# mesh.updateGeometry(disp)
mesh.plotDeformed(disp, 21.5)

def solve():

    # Initial stiffness matrix
    hbw = model.get_Matrix_0(mbuild, fint, disp, mesh)
    print("The half-band-width is", hbw)

    # K = mbuild.getBlock(range(20,25), range(20,25))
    # print(K)

    ndof = cons.dofCount()
    fdof = cons.get_fdof()

    # Initialize Data: Da = 0, r = fext - fint
    Da = np.zeros(ndof)
    da = np.zeros(ndof)

    r = fext - fint

    solver = Solver("numpy", cons)
    K = mbuild.getDenseMatrix()
    solver.solve(K, da, r, hbw)

    # K = mbuild.getMatrix()
    # plt.spy(K)
    # plt.show()

    # Update displacement vector
    Da[fdof] = Da[fdof] + da[fdof]
    print(Da)
    return Da
