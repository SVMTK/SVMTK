"""smooth a surface."""

import brainmesh as bm
from pathlib import Path

outdir = Path("results")
outdir.mkdir(exist_ok=True)

pig = bm.BrainSurface("data/pig.off")

for i in range(5):
    pig.smooth_laplacian(0.8)
    pig.smooth_laplacian(-0.805)        # Taubin smoothing

pig.smooth_taubin(5)

pig.save(str(outdir/"smoothed_pig.off"))
