# Import Standard Libraries
import matplotlib.pyplot as plt
import scipy as np
import time

# Import Local Libraries
from properties import Properties
from globalData import GlobalData
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

# Global Data
globdat = GlobalData()

# Mesh into Global Data
module = InputModule("input")
module.init(props, globdat)

# Global Data
module = InitModule()
module.init(props, globdat)

# Constraints
if dim == 2:
    idofs = globdat.mesh.getDofIndices(0, ["u", "v"])
    globdat.cons.addConstraints(idofs)
    idof = globdat.mesh.getDofIndex(2, "v")
    globdat.cons.addConstraint(idof)
    idof = globdat.mesh.getDofIndex(1, "u")
    globdat.cons.addConstraint(idof)

    # Loads
    idof = globdat.mesh.getDofIndex(1, "v")
    globdat.fext[idof] = 10

elif dim == 3:
    idofs = globdat.mesh.getDofIndices([0, 3], ["u", "v", "w"])
    globdat.cons.addConstraints(idofs)
    idofs = globdat.mesh.getDofIndices([1, 4], ["v", "w"])
    globdat.cons.addConstraints(idofs)
    idofs = globdat.mesh.getDofIndices([2, 5], ["u", "w"])
    globdat.cons.addConstraints(idofs)

    # Loads
    idofs = globdat.mesh.getDofIndices([2, 5], "v")
    globdat.fext[idofs] = 5

# Initial stiffness matrix
hbw = globdat.model.get_Matrix_0(globdat.mbuild, globdat.fint, globdat.disp, globdat.mesh)
print("The half-band-width is",hbw)

# Solve
solver = Solver("numpy", globdat.cons)
K = globdat.mbuild.getDenseMatrix()
solver.solve(K, globdat.disp, globdat.fext, hbw)

# K = mbuild.getMatrix()
# plt.spy(K)
# plt.show()

stop = time.time()
print("Elapsed time is ",stop-start)

# globdat.mesh.updateGeometry(disp)
globdat.mesh.plotDeformed(globdat.disp, 21.5)

def solve():

    # Initial stiffness matrix
    hbw = model.get_Matrix_0(mbuild, fint, disp, globdat.mesh)
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
