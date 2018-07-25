"""Remesh a surface."""

import brainmesh as bm
from pathlib import Path

outdir = Path("results")
outdir.mkdir(exist_ok=True)

pig = bm.BrainSurface("data/pig.off")
target_edge_length = 0.04
nb_iter = 3
protect_border = True

pig.isotropic_remeshing(
    target_edge_length,
    nb_iter,
    protect_border
)

pig.save(str(outdir/"remeshed_pig.off"))
