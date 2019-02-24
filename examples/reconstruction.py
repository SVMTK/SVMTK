"""Compute number of self-intersections."""

import SVMTK as  svm
from pathlib import Path

outdir = Path("results")
outdir.mkdir(exist_ok=True)

pig = svm.Surface("data/blobby.off")

self_intersecting = pig.self_intersections()
if self_intersecting:
    print("There are self-intersections")
    print(f"{pig.num_self_intersections()} pars of triangles intersect")
else:
    print("There are no self-intersecting")

pig.reconstruct_surface()

self_intersecting = pig.self_intersections()
if self_intersecting:
    print("There are self-intersections")
    print(f"{pig.num_self_intersections()} pars of triangles intersect")
else:
    print("There are no self-intersecting")

pig.save("foo.off")
