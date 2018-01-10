"""Script for converting between surface file formats."""
import os

from argparse import ArgumentParser

from .mesh_conversion_utils import (
    readAsc,
    readOff,
    readSTL,
    srf2off,
    srf2off_vec,
    srf2stl,
    write_xyz,
)


def brainmeshConvert(inpath: str, outpath: str) -> None:
    """Read `infile`, and save to `outfile`. Formats are inferred from suffixes.

    Arguments:
        infile: Input file name.
        outfile: output file name

    Supported input formats are:
     - .asc
     - .srf
     - .off

    Supported output formats are:
     - .stl
     - .off
     - .xyz
    """
    valid_insuffix = {".asc", ".srf", ".off", ".stl"}
    valid_outsuffix = {".stl", ".off", ".xyz"}

    insuffix = os.path.splitext(inpath)[-1]
    outsuffix = os.path.splitext(outpath)[-1]

    msg = "Invalid filetype, Expected {}".format(" or ".join(valid_insuffix))
    assert insuffix in valid_insuffix, msg
    msg = "Invalid filetype, Expected {}".format(" or ".join(valid_outsuffix))
    assert outsuffix in valid_outsuffix, msg

    if insuffix in {".asc", ".srf"}:
        datatuple = readAsc(inpath)
    elif insuffix == ".off":
        datatuple = readOff(inpath)
    elif insuffix == ".stl":
        datatuple = readSTL(inpath)

    if outsuffix == ".off":
        srf2off(datatuple.data, datatuple.num_vertices, datatuple.num_facets, outpath)
        # srf2off_vec(datatuple.data, datatuple.num_vertices, datatuple.num_facets, outpath)
    elif outsuffix == ".stl":
        srf2stl(datatuple.data, datatuple.num_vertices, datatuple.num_facets, outpath)
    elif outsuffix == ".xyz":
        write_xyz(datatuple.data, datatuple.num_vertices, outpath)


def createParser() -> ArgumentParser:
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
    parser = createParser()
    args = parser.parse_args()
    brainmeshConvert(args.input, args.output)


if __name__ == "__main__":
    main()
