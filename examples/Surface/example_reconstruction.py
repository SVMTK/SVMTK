"""Compute number of self-intersections."""

import SVMTK as  svm
from pathlib import Path

outdir = Path("results")
outdir.mkdir(exist_ok=True)

pig = svm.Surface("../Data/blobby.off")

self_intersecting = pig.num_self_intersections()
if self_intersecting>0:
    print("There are self-intersections")
    print("{} pars of triangles intersect".format(pig.num_self_intersections()))
else:
    print("There are no self-intersecting")

pig.reconstruct(20,0.1,0.1)

self_intersecting = pig.num_self_intersections()
if self_intersecting>0:
    print("There are self-intersections")
    print("{} pars of triangles intersect".format(pig.num_self_intersections()))
else:
    print("There are no self-intersecting")

pig.save("foo.off")
