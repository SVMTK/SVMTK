import numpy as np


def add_subdomain_tags(slice: "CGALSlice", cell_list, point_list) -> np.ndarray:
    """Add a subdomain tag to every cell.

    The tag is based on the smallest sice sontraint containting a cell vertex.

    NB! experimental!

    Arguments:
        slice: A CGALSlice instance.
        cell_list: list of cells.
        point_list: list of vertices.

    Returns:
        np.ndarray of mesh tags in the same order as cell_list.
    """
    cell_list = [
        max((slice.subdomain_map(x, y) for x, y, _ in point_list[cell])) for cell in cell_list
    ]
    return np.asarray(cell_list, dtype=np.uint64)
