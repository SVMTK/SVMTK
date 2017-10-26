"""Utilities for .asc or .srf and writing .off or .stl."""
import os
import sys 

from collections import namedtuple

from itertools import chain

import numpy as np


__STL_TEMPLATE = """
facet normal {:e} {:e} {:e}
    outer loop
        vertex {:e} {:e} {:e}
        vertex {:e} {:e} {:e}
        vertex {:e} {:e} {:e}
    endloop
endfacet\n
"""


DataTuple = namedtuple("DataTuple", ["num_facets", "num_vertices", "data"])


def facetNormals_vec(vertices: np.array) -> np.array:
    """Compute facet normals.
   
    Arguments:
        vertices: np.array of shape (N, 3, 3), where N is the number of facets

    Returns:
        (N, 3) unit outward facing facet normals.
    """
    v1 = vertices[:, 1] - vertices[:, 0]
    v2 = vertices[:, 2] - vertices[:, 0] 
    n = np.cross(v1, v2)
    n /= np.linalg.norm(n, axis=1)[:, None]
    return n


def srf2off(data: np.array, num_vertices: int, num_facets: int, outname: str) -> None:
    """Convert asc surface to off.

    Arguments:
        data: Array of vertices and facets of shape (N, 3).
        num_vertices: Number of vertices.
        num_facets: Number of facets.
    """
    assert data.shape[1] == 3                           # Number of vertices in a facet
    assert data.shape[0] == num_vertices + num_facets   # Header is consistent with data

    with open(outname, "w") as outfile:
        outfile.write("OFF\n")
        outfile.write("{Nv} {Nf} {Ne}".format(Nv=num_vertices, Nf=num_facets, Ne=0))

        # Write vertices
        for line in data[:num_vertices]:
            outfile.write("%e %e %e\n".format(*map(float, line)))

        # Write facet connectivity
        for line in data[num_vertices:]:
            outfile.write("3 %d %d %d\n".format(*map(int, line)))


def srf2off_vec(data: np.array, num_vertices: int, num_facets: int, outname: str) -> None:
    """Convert .asc surface to .off.

    Arguments:
        data: Array of vertices and facets of shape (N, 3).
        num_vertices: Number of vertices.
        num_facets: Number of facets.
    """
    assert data.shape[1] == 3                           # Numver of vertices in a facet
    assert data.shape[0] == num_vertices + num_facets   # Headr is consistent with data

    # Assume outfile does not exist, so we can append
    with open(outname, "ab") as outfile:
        # Write header
        outfile.write("OFF\n".encode())
        outfile.write("{Nv} {Nf} {Ne}\n".format(
            Nv=num_vertices, Nf=num_facets, Ne=0).encode()
        )
        # Use numpy to append vertex array
        np.savetxt(outfile, data[:num_vertices], fmt="%.15f")

        # Prepend with column of number of vertices
        facets = np.concatenate(
            (
                3*np.ones(num_facets, dtype=np.int16)[None].T,
                data[num_vertices:].astype(np.int16),
            ),
            axis=1
        )
        # Use numpy to append vertex connectivity array
        np.savetxt(outfile, facets, fmt="%d")


def srf2stl(data: np.array, num_vertices: int, num_facets: int, outname: str) -> None:
    """Convert asc surface to stl.

    Arguments:
        data: Array of vertices and facets of shape (N, 3).
        num_vertices: Number of vertices.
        num_facets: Number of facets.
    """
    name = os.path.splitext(os.path.basename(outname))[0]   # name for stl header
    with open(outname, "w") as outfile:
        outfile.write("solid {name}\n".format(name=name))

        # Use efancy indexing to get vertices. Probably super memory inefficient
        facets = data[data[num_vertices:].astype(np.int16)]

        # Compute normals
        normals = facetNormals_vec(facets)

        # Fill in template for each facet
        for i in range(num_facets):
            outfile.write(
                __STL_TEMPLATE.format(*chain(normals[i], facets[i].flatten()))
            )
        outfile.write("endsolid {name}".format(name=name))


def readAsc(filename: str) -> DataTuple:
    """Verify that `filename` has the correct format, and return the data.

    Arguments:
       filename: Name of inout file

    Returns:
        `DataTuple` with number of vertices, facets, the coordinates and the connectivity.
    """
    msg = "Cannot fine file {filename}".format(filename=filename)
    assert os.path.isfile(filename), msg

    data = np.loadtxt(filename, skiprows=2)
    msg = "Unrecognised file format. Expected 4 columns"
    assert data.shape[1] == 4, msg

    with open(filename, "r") as inputfile:
        line = inputfile.readline()     # Skip first line
        num_vertices, num_facets = tuple(map(int, inputfile.readline().split()))

        msg = "The number of vertices and facets is inconsistent with data size"
        assert data.shape[0] == num_vertices + num_facets, msg

    # Discard the column with number of edges (at least for .off)
    return DataTuple(num_facets=num_facets, num_vertices=num_vertices, data=data[:, :-1])


def readOff(filename: str) -> DataTuple:
    """Verify that `filename` has the correct format, and return the data.

    Arguments:
        filename: Name of input file (.off)

    Returns:
        `DataTuple` with number of vertices, facets, the coordinates and the connectivity.
    """
    msg = "Cannot fine file {filename}".format(filename=filename)
    assert os.path.isfile(filename), msg

    with open(filename, "r") as infile:
        assert "OFF" in infile.readline(), "Cannot find 'OFF' in header"
        num_vertices, num_facets, _ = tuple(map(int, infile.readline().split()))
        # TODO: Write a regex

        data_lines = list(map(lambda x: x.strip().split(), infile.readlines()))
        vertices = np.array(data_lines[:num_vertices], dtype=np.float32)
        facets = np.array(data_lines[num_vertices:], dtype=np.int16)[:, 1:]

    msg = "vertices.shape: {}, facets.shape: {}".format(vertices.shape, facets.shape)
    assert vertices.shape[1] == facets.shape[1], msg
    msg = "expected: {}, got {}".format(vertices.shape[0], num_vertices)
    assert vertices.shape[0] == num_vertices, msg
    msg = "expected: {}, got {}".format(facets.shape[0], num_facets)
    assert facets.shape[0] == num_facets, msg
    data = np.concatenate((vertices, facets))

    # data = np.loadtxt(filename, skiprows=2)

    msg = "Expected data.shape[0] == {0}, got {1}"
    msg.format(num_vertices + num_facets, data.shape[0])
    assert data.shape[0] == num_vertices + num_facets, msg

    return DataTuple(num_facets=num_facets, num_vertices=num_vertices, data=data)


def readSTL(filename: str) -> DataTuple:
    """Verify that `filename` has the correct format, and return the data.

    Arguments:
        filename: Name of input file (.stl)

    Returns:
        `DataTuple` with number of vertices, facets, the coordinates and the connectivity.
    """
    raise NotImplementedError
