"""Collapse long edges."""

import brainmesh as bm
from pathlib import Path

outdir = Path("results")
outdir.mkdir(exist_ok=True)

shark = bm.BrainSurface("data/mech-holes-shark.off")

print(f"Number of edges before: {shark.num_edges()}")
shark.collapse_edges(0.8)
print(f"Number of edges after: {shark.num_edges()}")

shark.save(str(outdir/"collapsed_shark.off"))
