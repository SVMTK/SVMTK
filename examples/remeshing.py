"""Remesh a surface."""

from svmtk import Surface
from pathlib import Path

outdir = Path("results")
outdir.mkdir(exist_ok=True)

pig = Surface("data/pig.off")
target_edge_length = 0.1
nb_iter = 3
protect_border = False

pig.isotropic_remeshing(
    target_edge_length,
    nb_iter,
    protect_border
)

pig.save(str(outdir / "remeshed_pig.off"))
