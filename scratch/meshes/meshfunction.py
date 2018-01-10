from dolfin import *
import numpy as np
from IPython import embed
import sys

mesh1 = UnitSquareMesh(100, 100)
mesh2 = UnitSquareMesh(150, 150)
mesh2.coordinates()[:] += 0.55

joined_cells = np.concatenate((mesh1.cells(), mesh2.cells() + mesh1.cells().max() + 1))
joined_cells = np.array(joined_cells, dtype=np.uintp)
joined_coordinates = np.concatenate((mesh1.coordinates(), mesh2.coordinates()))


mesh = Mesh("square.xml")
# mesh = Mesh()
# editor = MeshEditor()

# editor.open(mesh, 2, 2)  # top. and geom. dimension are both 2
# editor.init_vertices(joined_coordinates.shape[0])  # number of vertices
# editor.init_cells(joined_cells.shape[0])     # number of cells

# for i, c in enumerate(joined_coordinates):
#     editor.add_vertex(i, c)

# for i, c in enumerate(joined_cells):
#     editor.add_cell(i, c)

# editor.close()

# plot(mesh)

class Brain(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] < 0.7 + 1e-3 or x[1] < 0.7 + 1e-3


class Ventricle(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] > 0.55 - 1e-3 and x[1] > 0.55 - 1e-3 


class Interface(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 0.55) and x[1] > 0.55 - 1e-3 \
            or near(x[1], 0.55) and x[0] > 0.55 - 1e-3


brain = Brain()
ventricle = Ventricle()
interface = Interface()


mf = CellFunction("size_t", mesh)
mf.set_all(0)
brain.mark(mf, 1)       # MArk brain fisst, Overwrite with ventricle later
ventricle.mark(mf, 2)

ff = FacetFunction("size_t", mesh)
interface.mark(ff, 3)
plot(mf)

# plot(mf)
# input()
# # embed()
# sys.exit(1)

a0 = Constant(1.0)
a1 = Constant(0.01)
g_L  = Expression("- 10*exp(-pow(x[1] - 0.5, 2))", degree=2)
g_R = Constant(1.0)
f = Constant(1.0)

# Define function space and basis functions
V = FunctionSpace(mesh, "CG", 2)
u = TrialFunction(V)
v = TestFunction(V)

bcs = [DirichletBC(V, 5.0, DomainBoundary())]

ds = Measure("ds", domain=mesh, subdomain_data=ff)
dx = Measure("dx", domain=mesh, subdomain_data=mf)

F = inner(a0*grad(u), grad(v))*dx(1) + inner(a1*grad(u), grad(v))*dx(2) \
    - g_L*v*ds(1)

# Separate left and right hand sides of equation
a, L = lhs(F), rhs(F)

# Solve problem
u = Function(V)
solve(a == L, u, bcs)

# Plot solution and gradient
plot(u, title="u")
# plot(grad(u), title="Projected grad(u)")
input()
