"""smooth a surface."""

import SVMTK as svm
from pathlib import Path

outdir = Path("results")
outdir.mkdir(exist_ok=True)

pig = svm.Surface("../Data/pig.off")

for i in range(5):
    pig.smooth_laplacian(0.8,1)
    pig.smooth_laplacian(-0.805,1)        # Taubin smoothing

pig.smooth_taubin(5)

pig.save(str(outdir/"smoothed_pig.off"))
