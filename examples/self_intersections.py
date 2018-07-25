"""Compute number of self-intersections."""

import brainmesh as  bm
from pathlib import Path

outdir = Path("results")
outdir.mkdir(exist_ok=True)

pig = bm.BrainSurface("data/pig.off")
self_intersecting = pig.self_intersections()
if self_intersecting:
    print("There are self-intersections")
    print(f"{pig.num_self_intersections()} pars of triangles intersect")
