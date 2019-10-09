"""Fill the holes in a surface."""

import SVMTK as svm
from pathlib import Path

outdir = Path("../results")
outdir.mkdir(exist_ok=True)

mech_shark = svm.Surface("../Data/mech-holes-shark.off")
mech_shark.fill_holes()
mech_shark.save(str(outdir/"mech-filled-shark.off"))
