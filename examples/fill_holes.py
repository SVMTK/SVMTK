"""Fill the holes in a surface."""

import brainmesh as bm
from pathlib import Path

outdir = Path("results")
outdir.mkdir(exist_ok=True)

mech_shark = bm.BrainSurface("data/mech-holes-shark.off")
mech_shark.fill_holes()
mech_shark.save(str(outdir/"mech-filled-shark.off"))
