"""Compute number of self-intersections."""

import SVMTK as  svm
from pathlib import Path
outdir = Path("results")
outdir.mkdir(exist_ok=True)

sphere = svm.Surface("data/sphere_r6.off")

mesh = svm.Domain(sphere)

dx = 1.0
vector1 = [ (0,0,i*dx-3.) for i in range(6) ] 

vectors = [vector1]
mesh.add_polylines(vectors)

mesh.create_mesh()
mesh.save("der.mesh")

