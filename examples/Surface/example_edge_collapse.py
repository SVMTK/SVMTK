"""Collapse long edges."""

import SVMTK as svm
from pathlib import Path

outdir = Path("results")
outdir.mkdir(exist_ok=True)

shark = svm.Surface("../Data/mech-holes-shark.off")

print("Number of edges before: {}".format( shark.num_edges()) )
shark.collapse_edges(0.8)
print("Number of edges after: {}".format(shark.num_edges()))

shark.save(str(outdir/"collapsed_shark.off"))
