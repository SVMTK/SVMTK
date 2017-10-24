"""Create a scatterplot of an .xyz file.

Note that meshlab can also displpay .xyz files.
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

from mpl_toolkits.mplot3d import Axes3D

from argparse import ArgumentParser

from .mesh_conversion_utils import (
    read_xyz,
)


def scatterplot(data: np.array, figname: str, ndots: int=None) -> None:
    """Create a scatterplot of unstructured data.

    Arguments:
        data: The data.
        figname: The name (and format) of the output
        ndots (optional): The maximum number of points that will be displayed
    """
    sns.set()   # Make pretty
    assert data.shape[1] >= 3, "Expecting input to be 3D."""

    if ndots is not None:
        random_indices = np.random.randint(data.shape[0], size=int(ndots), dtype=np.int32)
        data = data[random_indices]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(data[:, 0], data[:, 1], data[:, 2])
    fig.savefig(figname)


def createParser() -> ArgumentParser:
    """Create argument parser with option (-i. --input) and optionally (-o, --output)."""
    parser = ArgumentParser(
        description="Save a scatterplot of the given .xyz file"
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
    parser.add_argument(
        "-n",
        "--ndots",
        help="Number of vertices to plot, drawn randomly",
        required=False
    )
    return parser


def main() -> None:
    """Execute the parser and call scatterplot."""
    parser = createParser()
    args = parser.parse_args()

    msg = "Expected .xyz input file."
    assert args.input.endswith(".xyz"), msg

    data = read_xyz(args.input)
    scatterplot(data, args.output, args.ndots)
