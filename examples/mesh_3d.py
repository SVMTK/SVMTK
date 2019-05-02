"""Create a 3D mesh"""

from svmtk import Surface, Domain
from pathlib import Path

outdir = Path("results")
outdir.mkdir(exist_ok=True)

blob = Surface("data/elephant.off")
print("num_vertices: ", blob.num_vertices())

mc = Domain(blob)
mc.create_mesh(16)

mc.save("results/3d_mesh")
