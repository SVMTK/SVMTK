"""Collapse long edges."""

import SVMTK as svm
from pathlib import Path

outdir = Path("results")
outdir.mkdir(exist_ok=True)

shark = svm.Surface("data/mech-holes-shark.off")

print(f"Number of edges before: {shark.num_edges()}")
shark.collapse_edges(0.8)
print(f"Number of edges after: {shark.num_edges()}")

shark.save(str(outdir/"collapsed_shark.off"))
