"""Convert between mesh formats.

Can read asc and srf files and falls back to meshio for everything else.
"""

from argparse import ArgumentParser
from pathlib import Path
import meshio

import numpy as np


def readAsc(filename: Path) -> meshio.Mesh:
    """Verify that `filename` has the correct format, and return the data.

    Arguments:
       filename: Name of inout file

    Returns:
        `DataTuple` with number of vertices, facets, the coordinates and the connectivity.
    """
    cols = (0, 1, 2)
    data = np.loadtxt(filename, skiprows=2, usecols=cols)

    msg = "Unrecognised file format. Expected 4 columns"
    assert data.shape[1] == len(cols), msg

    with open(filename, "r") as inputfile:
        line = inputfile.readline()     # Skip first line
        num_vertices, num_facets = tuple(map(int, inputfile.readline().split()))

        msg = "The number of vertices and facets is inconsistent with data size"
        assert data.shape[0] == num_vertices + num_facets, msg

    vertices = data[:num_vertices]
    facets = data[num_vertices:].astype(np.int64)
    datatuple = meshio.Mesh(
        points=vertices,
        cells={"triangle": facets}
    )
    return datatuple


def brainmeshConvert(input_name: str, output_name: str) -> None:
    """Read `infile`, and save to `outfile`. Formats are inferred from suffixes.

    Arguments:
        input_name: Input file name.
        output_name: output file name
    """
    inpath = Path(input_name)
    outpath = Path(output_name)

    msg = "Cannot find file: {input_name}".format(input_name=input_name)
    assert inpath.exists(), msg

    if inpath.suffix in {".asc", ".srf"}:
        mesh = readAsc(inpath)
    else:
        mesh = meshio.read(str(inpath))
    meshio.write(str(outpath), mesh)


def create_parser() -> ArgumentParser:
    """Create argument parser with option (-i. --input) and optionally (-o, --output)."""
    parser = ArgumentParser(
        description="Convert surface file from .asc to .off or .stl"
    )
    parser.add_argument(
        "-i",
        "--input",
        help="path to input file",
        required=True
    )
    parser.add_argument(
        "-o",
        "--output",
        help="path to output file",
        required=True
    )
    return parser


def main() -> None:
    """brainmesh-convert entry point."""
    parser = create_parser()
    args = parser.parse_args()
    brainmeshConvert(args.input, args.output)


if __name__ == "__main__":
    main()
