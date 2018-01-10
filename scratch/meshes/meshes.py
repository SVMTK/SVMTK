from dolfin import *
import numpy as np


mesh1 = UnitSquareMesh(10, 10)
mesh2 = UnitSquareMesh(15, 15)
mesh2.coordinates()[:] += 0.5

joined_cells = np.concatenate((mesh1.cells(), mesh2.cells() + mesh1.cells().max() + 1))
joined_cells = np.array(joined_cells, dtype=np.uintp)
joined_coordinates = np.concatenate((mesh1.coordinates(), mesh2.coordinates()))


mesh = Mesh()
editor = MeshEditor()

editor.open(mesh, 2, 2)  # top. and geom. dimension are both 2
editor.init_vertices(joined_coordinates.shape[0])  # number of vertices
editor.init_cells(joined_cells.shape[0])     # number of cells

for i, c in enumerate(joined_coordinates):
    editor.add_vertex(i, c)

for i, c in enumerate(joined_cells):
    editor.add_cell(i, c)

editor.close()
plot(mesh)
input()

V = FunctionSpace(mesh, "CG", 1)

u = TrialFunction(V)
v = TestFunction(V)

m = u*v*dx

f = Expression("x[0] + x[1]", degree=1)*v*dx

U_ = Function(V)

M = assemble(m)
b = assemble(f)

solve(M, U_.vector(), b, "cg", "jacobi")

plot(U_)
input()
