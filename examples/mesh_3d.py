"""Create a 3D mesh"""

import brainmesh as bm
from pathlib import Path

outdir = Path("results")
outdir.mkdir(exist_ok=True)

blob = bm.BrainSurface("data/elephant.off")
print("num_vertices: ", blob.num_vertices())

mc = bm.BrainMesh(blob)
mc.create_mesh()
# mc.refine_mesh()

mc.save("results/3d_mesh")
