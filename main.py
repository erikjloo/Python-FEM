# Import Standard Libraries
import scipy as np

# Import Local Libraries
from properties import Properties
from constraints import Constraints
from algebra import MatrixBuilder
from modules import InputModule
from solvers import Solver
from models import ModelFactory
from mesh import Mesh

# Initialization
file = "Examples/semicircle.pro"
props = Properties()
props.parseFile(file)
props.print()

# Mesh
module = InputModule("input")
mesh = module.init(props)

# Model
model = ModelFactory("model", props, mesh)

# multistep(model, mesh, nsteps=1, nrkey=None)

# Global Data
ndof = mesh.dofCount()
cons = Constraints(ndof)
mbuild = MatrixBuilder(ndof)
fint = np.zeros(ndof)
fext = np.zeros(ndof)
disp = np.zeros(ndof)

# Constraints
idofs = mesh.getDofIndices([0, 3], ["u", "v", "w"])
cons.addConstraints(idofs)
idofs = mesh.getDofIndices([1, 4], ["v", "w"])
cons.addConstraints(idofs)
idofs = mesh.getDofIndices([2, 5], ["u", "w"])
cons.addConstraints(idofs)

# Loads
idofs = mesh.getDofIndices([2, 5], ["v"])
fext[idofs] = 5

# Initial stiffness matrix
model.get_Matrix_0(mbuild, fint, disp, mesh)
K = mbuild.getBlock(range(20,25), range(20,25))

ndof = cons.dofCount()
fdof = cons.get_fdof()

# Initialize Data: Da = 0, r = fext - fint
Da = np.zeros(ndof)
da = np.zeros(ndof)

r = fext - fint
max_hbw = 0

solver = Solver("python", cons)
K = mbuild.getDenseMatrix()
solver.solve(K, da, r)
print(da)
# Update displacement vector
Da[fdof] = Da[fdof] + da[fdof]
# mesh.updateGeometry(d+Da)

mesh.plotDeformed(Da, 37)

# Modules works as a series of commands

def execute(module):
    module.init()
    module.run()
    module.shutdown()

# if __name__ == "__main__":
#     execute(module)
