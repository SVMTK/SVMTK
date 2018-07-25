"""Convert a surface to a triangular mesh."""

import brainmesh as bm
from pathlib import Path

outdir = Path("results")
outdir.mkdir(exist_ok=True)

p = bm.BrainSurface("data/P.off")
p.triangulate_faces()
p.save(str(outdir/"triangularP.off"))
