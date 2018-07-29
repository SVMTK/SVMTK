"""Compute boolean operations between blobby and eight."""

import brainmesh as bm
from pathlib import Path


def print_stats(surf: bm.BrainSurface) -> None:
    """Print num_facet, num_edges and num vertices."""
    print(f"num_faces: {surf.num_faces()}")
    print(f"num_edges: {surf.num_edges()}")
    print(f"num_vertices: {surf.num_vertices()}")


class BooleanOperation:
    """Class for performing boolean operation on surffaces.

    Operations include intersection, union and difference.
    """

    def __init__(self, input_name1: str, input_name2: str) -> None:
        """Load surfaces from filepaths."""
        self.surf1 = bm.BrainSurface(input_name1)
        self.surf2 = bm.BrainSurface(input_name2)

    def intersection(self) -> None:
        """Copute the intersection between surface1 and surface2."""
        self.surf1.intersection(self.surf2)

    def union(self) -> None:
        """Copute the union between surface1 and surface2."""
        self.surf1.union(self.surf2)

    def difference(self) -> None:
        """Copute the difference between surface1 and surface2."""
        self.surf1.difference(self.surf2)

    def save(self, output_name: str) -> None:
        """Save surface1 as `output_name`."""
        self.surf1.save(str(output_name))       # NB! Must convert to string for pybind11


if __name__ == "__main__":
    outdir = Path("results")
    outdir.mkdir(exist_ok=True)

    bop = BooleanOperation(
        "data/blobby.off",
        "data/eight.off"
    )

    print_stats(bop.surf1)
    bop.difference()
    bop.save(outdir/"blobby_difference_eight.off")

    # Operations are performed in-place on surface1
    bop = BooleanOperation(
        "data/blobby.off",
        "data/eight.off"
    )

    bop.union()
    bop.save(outdir/"blobby_union_eight.off")

    # 
    bop = BooleanOperation(
        "data/blobby.off",
        "data/eight.off"
    )

    bop.intersection()
    bop.save(outdir/"blobby_intersection_eight.off")
