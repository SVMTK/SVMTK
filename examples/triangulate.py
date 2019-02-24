"""Convert a surface to a triangular mesh."""

import SVMTK as svm
from pathlib import Path

outdir = Path("results")
outdir.mkdir(exist_ok=True)

p = svm.Surface("data/P.off")
p.triangulate_faces()
p.save(str(outdir/"triangularP.off"))
