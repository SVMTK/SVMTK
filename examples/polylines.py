"""Compute number of self-intersections."""

import brainmesh as  bm
from pathlib import Path
outdir = Path("results")
outdir.mkdir(exist_ok=True)

sphere = bm.BrainSurface("data/sphere_r6.off")

mesh = bm.BrainMesh(sphere)

dx = 1.0
vector1 = [ (0,0,i*dx-3.) for i in range(6) ] 

vectors = [vector1]
mesh.add_polylines(vectors)

mesh.create_mesh()
mesh.save("der.mesh")

