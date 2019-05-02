"""Convert a surface to a triangular mesh."""

from svmtk import Surface
from pathlib import Path

outdir = Path("results")
outdir.mkdir(exist_ok=True)

p = Surface("data/P.off")
p.triangulate_faces()
p.save(str(outdir/"triangularP.off"))
