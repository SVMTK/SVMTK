/*
  This file contains docstrings for use in the Python bindings.
  Do not edit! They were automatically extracted by pybind11_mkdoc.
 */

#define __EXPAND(x)                                      x
#define __COUNT(_1, _2, _3, _4, _5, _6, _7, COUNT, ...)  COUNT
#define __VA_SIZE(...)                                   __EXPAND(__COUNT(__VA_ARGS__, 7, 6, 5, 4, 3, 2, 1))
#define __CAT1(a, b)                                     a ## b
#define __CAT2(a, b)                                     __CAT1(a, b)
#define __DOC1(n1)                                       __doc_##n1
#define __DOC2(n1, n2)                                   __doc_##n1##_##n2
#define __DOC3(n1, n2, n3)                               __doc_##n1##_##n2##_##n3
#define __DOC4(n1, n2, n3, n4)                           __doc_##n1##_##n2##_##n3##_##n4
#define __DOC5(n1, n2, n3, n4, n5)                       __doc_##n1##_##n2##_##n3##_##n4##_##n5
#define __DOC6(n1, n2, n3, n4, n5, n6)                   __doc_##n1##_##n2##_##n3##_##n4##_##n5##_##n6
#define __DOC7(n1, n2, n3, n4, n5, n6, n7)               __doc_##n1##_##n2##_##n3##_##n4##_##n5##_##n6##_##n7
#define DOC(...)                                         __EXPAND(__EXPAND(__CAT2(__DOC, __VA_SIZE(__VA_ARGS__)))(__VA_ARGS__))

#if defined(__GNUG__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#endif


static const char *__doc_DefaultMap =
R"doc(The defualt method to set subdomains in the mesh. Uses
bitstring to integer conversion to set subdomain tag.)doc";


static const char *__doc_DefaultMap_index =
R"doc(Maps bistring to an integer using binary conversion to integer

:param bits: bitstring of the form boost::dynamic_bitset.

:Returns: Float of the bit string.

)doc";

static const char *__doc_DefaultMap_make_interfaces =
R"doc(Gives each interfaces an unique tag based on presence in mesh and from the highest cell tag. Ex:(Cell tags 1 2 -> Interface tags 3 4 5)

:param interfaces: List of tuples integers, representing the interface between two cell tags. 

:Returns: Dictonary with tuple of two integer as key and an unique integer as value.

)doc";


static const char *__doc_Domain =
R"doc(The SVMTK Domain class is used to create and tetrahedra mesh by utilizing `CGAL <https://www.cgal.org/>`_. 

The SVMTK Domain class was implemented with the `Exact predicates inexact constructions kernel <https://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Exact__predicates__inexact__constructions__kernel.html>`_.  The construction of tetrahedra mesh in CGAL require `CGAL Domain <https://doc.cgal.org/latest/Mesh_3/group__PkgMesh3Domains.html>`_. 

Input surface are used to create `Polyhedral_mesh_domain_with_features_3 <https://doc.cgal.org/latest/Mesh_3/classCGAL_1_1Polyhedral__mesh__domain__with__features__3.html>`_. which are function wrapped to create a `Labeled_Mesh_Domain <https://doc.cgal.org/latest/Mesh_3/classCGAL_1_1Labeled__mesh__domain__3.html>`_. This allows the construction of meshes to create cell tags so that we may speicify different subdomains. The Labeled mesh domain properties are then used to declare `Mesh_domain_with_polyline_features_3 <https://doc.cgal.org/latest/Mesh_3/classCGAL_1_1Mesh__domain__with__polyline__features__3.html>`_.

`Mesh_Criteria <https://doc.cgal.org/latest/Mesh_3/classMeshCriteria__3.html>`_  defines construction process. The tetrahedra mesh is stored in `Mesh_complex_3_in_triangulation_3 <https://doc.cgal.org/latest/Mesh_3/classMeshComplex__3InTriangulation__3.html>`_ structure, where the data follows `Mesh_triangulation_3 <https://doc.cgal.org/latest/Triangulation_3/index.html>`_ triganulation data structure.

Notable CGAL declarations: |
- `Exact predicates inexact constructions kernel <https://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Exact__predicates__inexact__constructions__kernel.html>`_.
- `Labeled_Mesh_Domain <https://doc.cgal.org/latest/Mesh_3/classCGAL_1_1Labeled__mesh__domain__3.html>`_.
- `Mesh_domain_with_polyline_features_3 <https://doc.cgal.org/latest/Mesh_3/classCGAL_1_1Mesh__domain__with__polyline__features__3.html>`_.
- `Mesh_complex_3_in_triangulation_3 <https://doc.cgal.org/latest/Mesh_3/classMeshComplex__3InTriangulation__3.html>`_.
- `Mesh_Criteria <https://doc.cgal.org/latest/Mesh_3/classMeshCriteria__3.html>`_.


Attributes:
- **c3t3** - `Mesh_complex_3_in_triangulation_3 <https://doc.cgal.org/latest/Mesh_3/classMeshComplex__3InTriangulation__3.html>`_ object.
- **borders** - list of sequential points that does not overlapp in a non-conforming way.
- **features** - list of sequential points.

)doc";

static const char *__doc_Domain_Domain =
R"doc(Constructor for meshing a single surface

:param surface: :class:`Surface` object.
:param error_bound: the error bound of the surface representation.

)doc";

static const char *__doc_Domain_Domain_2 =
R"doc(Constructor for meshing multiple surfaces

:param surfaces: List of :class:`Surface` objects.
:param error_bound: the error bound of the surface representation.

)doc";

static const char *__doc_Domain_Domain_3 =
R"doc(Constructor for meshing multiple surfaces

:param surfaces: List of :class:`Surface` objects.
:param map: :class:`SubDomainMap` object, used to set subdomain and boundary tags.
:param error_bound: the error bound of the surface representation.

)doc";


static const char *__doc_Domain_Domain_4 =
R"doc(Constructor for loading a volumetric mesh from file.

Note: Supports only .mesh format with integer subdomain tag, and will rebind facet tags.

:param filename: the filename of the volumetric mesh. 
:param error_bound: the error bound of the surface representation.
:param

)doc";

static const char *__doc_Domain_add_border =
R"doc(Adds a polyline to the Domain attribute borders.

:param polyline: List of sequential points that represent an edges in the mesh. 

)doc";

static const char *__doc_Domain_add_feature =
R"doc(Adds a polyline to the Domain attribute features.

:param polyline: List of sequential points that represent an edges in the mesh. 

)doc";

static const char *__doc_Domain_clear_borders = R"doc(Clear borders.)doc";

static const char *__doc_Domain_clear_features = R"doc(Clear features.)doc";

static const char *__doc_Domain_get_subdomains =
R"doc(Returns a set of integer that represents the cell tags in the mesh.

:Returns: List of integers that represents all the cell tags in the mesh.

)doc";


static const char *__doc_Domain_get_facets =
R"doc(Returns all triangles as a numpy array (N,3) with vertex indices. The first vertex index starts at 1.  

:Returns: a numpy array of facets as vertex indices.

)doc";

static const char *__doc_Domain_get_facet_tags =
R"doc(Returns all facet tags as a numpy array, with the option of excluding unmarked facets.  

:param exclude_unmarked: boolean option to either exclude facets with tag equal to zero.

:Returns: a numpy array of facet tags.

)doc";

static const char *__doc_Domain_get_cells =
R"doc(Returns all tetrahedrons as a numpy array (N,4) with vertex indices. The first vertex index starts at 1.  

:Returns: a numpy array of cells as vertex indices.

)doc";

static const char *__doc_Domain_get_cell_tags =
R"doc(Returns all tetrahedron tags as a numpy array.   

:Returns: a numpy array with tetrahedron tags.

)doc";


static const char *__doc_Domain_get_features =
R"doc(Returns the added (internal) features to the :class:`Domain`.

:Returns: an array with :class:`Point_3`

)doc";

static const char *__doc_Domain_get_borders =
R"doc(Returns the added borders to the :class:`Domain`.

:Returns: an array with :class:`Point_3`

)doc";


static const char *__doc_Domain_refine =
R"doc(Refine the mesh to a target mesh_resolution.

:param mesh_resolution: determines the ratio between the minimum bounding radius and the upper bound of the cell circumradi.

)doc";


static const char *__doc_Domain_get_points =
R"doc(Returns the Cartesian coordinates for all vertices as a numpy array (N,3) with vertex locations.   

:Returns: a numpy array of all vertex coordinates.

)doc";



static const char *__doc_Domain_tetrahedral_remeshing =
R"doc( Tetrahedral remeshing of the tetrahedral mesh.

See: `Tetrahedral remeshing <https://doc.cgal.org/latest/Tetrahedral_remeshing/index.html>`_.

:param edge_length: is edge length after remeshing.
:param nb_iter: the number of iterations for the sequence of atomic operations performed. 
:param protect_borders: If true, constraint edges cannot be modified at all during the remeshing process. 

:Returns: None 

)doc";

static const char *__doc_Domain_add_sharp_border_edges =
R"doc(Adds sharp border edges from a SVMTK Surface object.

Checks polyhedron for sharp edges, if these edges do not conflict/intersect previously stored edges, then the edges are stored in member variable borders. The use of 1D features in combination with ill-posed meshing parameters can cause segmentation fault (crash).

:param surface: :class:`Surface` object.
:param threshold: the threshold angle in degree (0-90). Angles between connected edges above the threshold angle is considered sharp.

)doc";

static const char *__doc_Domain_add_sharp_border_edges_2 =
R"doc(Adds sharp border edges from a SVMTK Surface object in a given plane. 

:param surface: :class:`Surface` object.
:param plane: excludes all edges that is not within a given error of the plane.
:param threshold: the threshold angle in degree (0-90). Angles between connected edges above the threshold angle is considered sharp.

)doc";

static const char *__doc_Domain_add_sharp_border_edges_3 =
R"doc(Adds sharp border edges from a SVMTK Surface object on the surface of another :class:`Surface`. 

:param surface: :class:`Surface` object.
:param clip: :class:`Surface` object used to excludes all edges that is not within a given error of the clip surface.
:param threshold: the threshold angle in degree (0-90). Angles between connected edges above the threshold angle is considered sharp.

)doc";

static const char *__doc_Domain_boundary_segmentations =
R"doc(Segments the boundary of a specified subdomain tag. Calls :func:`Surface.surface_segmentation` on the boundary surfaces specific by input, and updates the boundary facets with the segmentation tags.

:param subdmain_tag: An integer representing a subdomain in the mesh. 
:param angle_in_degree: The threshold angle in degree used to detect sharp edges (0-90).

)doc";

static const char *__doc_Domain_boundary_segmentations_2 =
R"doc(

:param subdmain_tag: An integer representing a subdomain in the mesh. 
:param angle_in_degree: The threshold angle in degree used to detect sharp edges (0-90).

)doc";

static const char *__doc_Domain_boundary_segmentations_3 =
R"doc(Segments the boundary for each subdomain in the stored mesh.

:param angle_in_degree: The threshold angle in degree used to detect sharp edges (0-90).

)doc";

static const char *__doc_Domain_check_mesh_connections =
R"doc(Checks the connections in the mesh for bad vertices and bad edges, and prints out the number of bad vertices and bad edges.

A bad vertex is a boundary vertex that is shared by two non-adjacent cells A bad edge is a boundary edge that is shared by more than 2 facets.

In order to remove bad edges and bad edges, the surfaces must be free of self intersection and have mesh resolution equal to the highest mesh resolution of the surfaces. The separation of narrow gaps and the thickening of thin mesh segments may also be required.  

)doc";

static const char *__doc_Domain_create_mesh =
R"doc(Creates the mesh stored in the class attribute c3t3 using `CGAL mesh criteria parameters <https://doc.cgal.org/latest/Mesh_3/classCGAL_1_1Mesh__criteria__3.html>`_.

:param edge_size: Sets the upper bound for the lengths of curve edges.
:param cell_size: Sets the upper-bound for the circumradii of the mesh tetrahedra.
:param facet_size: Sets the upper-bound or for the radii of the surface Delaunay balls.
:param facet_angle: Sets the lower-bound for the angles (in degrees) of the surface mesh facets.
:param facet_distance: Sets the upper-bound for the distance between the facet circumcenter and the center of its surface Delaunay ball.
:param cell_radius_edge_ratio: Sets the upper-bound for the radius-edge ratio of the mesh tetrahedra.

)doc";


static const char *__doc_Domain_create_mesh_2 =
R"doc(Creates the mesh stored in the class attribute c3t3. The mesh criteria is set based on mesh_resolution and the minimum bounding radius of the mesh.

:param mesh_resolution: Sets the mesh criteria parameters: cell_size, facet_size, edge_size and facet_distance, by dividing the minimum bounding radius of the mesh with the mesh resolution. 

)doc";


static const char *__doc_Domain_create_mesh_3 =
R"doc(Creates the mesh stored in the class attribute c3t3. The mesh resolution is set automatic to the highest mesh resolution of the surface(s).

)doc";

static const char *__doc_Domain_init_triangulation =
R"doc(Inserts all points from a class:`Surface` object into the volumetric mesh before any triangulation.

Note: 
     Experimental, may result in failure.
     
:param surf: a class:`Surface` object:     
      

)doc";

static const char *__doc_Domain_init_triangulation_2 =
R"doc(Inserts all points from a list of  class:`Surface` objects into the volumetric mesh before any triangulation.

Note: 
     Experimental, may result in failure.
     
:param surf: a list of class:`Surface` object:     
      

)doc";


static const char *__doc_Domain_dihedral_angles =
R"doc(Computes the dihedral angle for all tetrahedron cells in complex (c3t3).

:Returns: List of dihedral angles in degree.

)doc";

static const char *__doc_Domain_dihedral_angles_min_max =
R"doc(Returns the maximum and minimum dihedral angles of cells in mesh.

:Returns: Tuple of floats with first as the minimum and second as maximum.

)doc";

static const char *__doc_Domain_exude =
R"doc(CGAL function for excude optimization of the constructed mesh.

See also:
`excude_optimize_mesh <https://doc.cgal.org/latest/Mesh_3/group__PkgMesh3Functions.html>`_.

:param time_limit: Sets, in seconds, a CPU time limit after which the optimization process is stopped.
:sliver_bound: Sets a targeted lower-bound on dihedral angles of mesh cells.

)doc";

static const char *__doc_Domain_get_collision_distances =
R"doc( Experimental: This function computes the collision distance in the negative normal for each facet on the interface for a subdomain. If
the collision occurs with a specified interface, then the collision distance is set as negative. The collisions distance is stored as 
triangle data, and can be written to file.

The purpose of this function is to approximate to the hydrolic resistance for a subdomain, then remove the 
subdomain, and simulate flow on the interface. 

:param subdomain_tag: An integer representing a subdomain in the mesh.  

:int boundary_tag=0:  An integer, which in combination wit the subdomain_tag gives an iterface in the mesh.


)doc";

static const char *__doc_Domain_get_collision_spheres =
R"doc(Experimental: This function computes the collision spheres for each facet on the interface for a subdomain. 
The collisions distance is stored as triangle data, and can be written to file.

The purpose of this function is to approximate to the hydrolic resistance for a subdomain, then remove the 
subdomain, and simulate flow on the interface. 

:param subdomain_tag: An integer representing a subdomain in the mesh.  


)doc";


static const char *__doc_Domain_write_facet_data =
R"doc(Writes stored triangle data to file.

:param filename: name of file. 


)doc";




static const char *__doc_Domain_get_interface =
R"doc(Returns the interface between two subdomain as a :class:`Surface` object.

:param interface: a pair of integers tags. 

:Returns: a :class:`Surface` object.
    
)doc";

static const char *__doc_Domain_get_boundaries =
R"doc(Iterates over all surface boundaries of subdomains and stores and returns it as a list of Surface objects.

:Returns: List of :class:`Surface` objects.
    
)doc";

static const char *__doc_Domain_get_boundary =
R"doc(Returns the boundary of a subdomain tag stored in a Surface object.

:param tag: An integer representing a subdomain in the mesh.  

:Returns: :class:`Surface` object.

)doc";

static const char *__doc_Domain_get_curve_tags =
R"doc(Returns the set of the curve/edge tags in the triangulation.

:Returns: set of integers representing all the curve tags.

)doc";

static const char *__doc_Domain_get_patches =
R"doc(Returns a set of integer pairs that represents the surface facet tags in the mesh.

:Returns: List of integer tuples that represents all the facet tags in the mesh.


)doc";

static const char *__doc_Domain_lloyd =
R"doc(CGAL function for lloyd optimization of the constructed mesh.

See also: `lloyd_optimize_mesh <https://doc.cgal.org/latest/Mesh_3/group__PkgMesh3Functions.html>`_. 

Note: 
Do not use lloyd after excude optimazation, it will cause failure.

:param time_limit: Sets , in seconds, a CPU time limit after which the optimization process is stopped.
:param max_iter: Sets a limit on the number of iterations.
:param convergence: the displacement of any vertex is less than a given percentage of the length of the shortest edge incident to that vertex.
:param freeze_bound: vertex that has a displacement less than a given percentage of the length (the of its shortest incident edge, is frozen (i.e. is not relocated).
:param do_freeze: - completes the freeze_bound parameter.

)doc";

static const char *__doc_Domain_num_cells =
R"doc(Returns number of tetrahedron cells in stored volume mesh.

:Returns: The number of tetrahedron cells in stored volume mesh.

)doc";

static const char *__doc_Domain_num_curves =
R"doc(Returns number of 1D polylines in stored volume mesh.

:Returns: The number of 1D polylines in stored volume mesh.

)doc";

static const char *__doc_Domain_num_facets =
R"doc(Returns number of facets, triangles, in stored volume mesh.

:Returns: The number of facets, triangles, in stored volume mesh.

)doc";

static const char *__doc_Domain_num_patches =
R"doc(Returns number of interface tags in stored volume mesh.

:Returns: The number of interface tags in stored volume mesh.

)doc";

static const char *__doc_Domain_num_subdomains =
R"doc(Returns number of subdomains in stored volume mesh.

:Returns: The number of subdomains in stored volume mesh.

)doc";

static const char *__doc_Domain_num_surfaces =
R"doc(Returns number of surfaces that was added in the constructor.

:Returns: The number of surfaces that was added in the constructor.

)doc";

static const char *__doc_Domain_num_vertices =
R"doc(Returns number of vertices in stored volume mesh.

:Returns: The number of vertices in stored volume mesh.

)doc";

static const char *__doc_Domain_odt =
R"doc(CGAL function for odt optimization of the constructed mesh.

See also: `odt_optimize_mesh <https://doc.cgal.org/latest/Mesh_3/group__PkgMesh3Functions.html>`_.

:param time_limit: Sets, in seconds, a CPU time limit after which the optimization process is stopped.
:param max_iter: Sets a limit on the number of performed iterations.
:param convergence: the displacement of any vertex is less than a given percentage of the length of the shortest edge incident to that vertex.
:param freeze_bound: vertex that has a displacement less than a given percentage of the length (the of its shortest incident edge, is frozen (i.e. is not relocated).
:paramdo_freeze: completes the freeze_bound parameter.

)doc";

static const char *__doc_Domain_perturb =
R"doc(CGAL function for perturb optimization of the constructed mesh.

See also: `perturb_mesh <https://doc.cgal.org/latest/Mesh_3/group__PkgMesh3Functions.html>`_.

:param time_limit: Sets, in seconds, a CPU time limit after which the optimization process is stopped.
:param sliver_bound: Sets targeted lower-bound on dihedral angles of mesh cells.

)doc";

static const char *__doc_Domain_radius_ratios =
R"doc(Computes the radius ratio for all tetrahedron cells in complex (c3t3).

:Returns: List of radius ratios.

)doc";

static const char *__doc_Domain_radius_ratios_min_max =
R"doc(Returns the maximum and minimum radius ratio of cells in mesh.

:Returns: Tuple of floats with first as the minimum and second as maximum.

)doc";

static const char *__doc_Domain_remove_subdomain =
R"doc(Removes all cells in the mesh with a specified integer tag, but perserves the interface tags as if no cells were removed.

:param tag: Removes cells with this specific integer tag. 
 
)doc";

static const char *__doc_Domain_remove_subdomain_2 =
R"doc(Removes all cells in the mesh with tags in a list, but perserves the interface tags as if no cells were removed.

:param tags: List of cell tag to be removed.

)doc";

static const char *__doc_Domain_save =
R"doc(Writes the mesh stored in the class attribute c3t3 to file.

The interface tags are loaded from :class:`SubdomainMap` added in the constructor. If there are no interfaces in :class:`SubDomainMap`, then default interfaces are selected.

:param filename: The path to the output file. 
:param save_1Dfeatures: Option to save the edges with tags.

)doc";


static const char *__doc_Plane3 =
R"doc(Wrapper for `CGAL Plane_3 class <https://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Plane__3.html>`_, with plane equation defined as :math:`h : ax+by+cz+d = 0.`  
        
)doc";

static const char *__doc_Plane3_Plane3 =
R"doc(Creates a SVMTK Plane_3 object by setting the coefficients of the plane equation.

:param a: Sets the first coefficient of the plane equation. 
:param b: Sets the second coefficient of the plane equation.
:param c: Sets the third coefficient of the plane equation.
:param d: Sets the fourth coefficient of the plane equation.        

)doc";

static const char *__doc_Plane3_Plane3_2 =
R"doc(Creates a SVMTK Plane_3 object by providing a point on the plane and the plane normal vector.

:param point: :class:`Point_3` object of a point on the plane.
:param vector: :class:`Vector_3` of the plane normal. 
        
)doc";

static const char *__doc_Point2 =
R"doc(Wrapper for `CGAL Point_2 class <https://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Point__2.html>`_.)doc";


static const char *__doc_Point2_Point2 =
R"doc(Creates a SVMTK Point_2 object, i.e. a point in xy plane.

:param x: Sets the x-coordinate of the point.
:param y: sets the y-coordinate of the point.
        
)doc";


static const char *__doc_Point3 =
R"doc(Wrapper for `CGAL Point_3 class <https://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Point__3.html>`_.)doc";

static const char *__doc_Point3_Point3 =
R"doc(Creates a SVMTK Point_3 object.

:param x: Sets the x coordinate of the point.
:param y: Sets the y coordinate of the point.
:param z: Sets the z coordinate of the point.        

)doc";


static const char *__doc_Slice =
R"doc(The SVMTK Slice class stores and manipulate triangulated surfaces in a plane, i.e. the third coordinate is neglected. The Slice class uses the `Exact predicates inexact constructions kernel <https://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Exact__predicates__inexact__constructions__kernel.html>`_. The meshing is done by triangulation of constrainted with edges by using CGAL class `Constrained_triangulation_plus_2 <https://doc.cgal.org/latest/Triangulation_2/classCGAL_1_1Constrained__triangulation__plus__2.html>`_ . (REWRITE) 
The Slice does not handle cavities, but cavities can be assigned with adding surfaces with :func:`add_surface_domains` and optional SVMTK :class:`SubdomainMap` object.

Notable CGAL declarations: 
- `Constrained_Delaunay_triangulation_2 <https://doc.cgal.org/latest/Triangulation_2/classCGAL_1_1Constrained__Delaunay__triangulation__2.html>`_.
- `Constrained_triangulation_plus_2 <https://doc.cgal.org/latest/Triangulation_2/classCGAL_1_1Constrained__triangulation__plus__2.html>`_
- `Delaunay_mesh_size_criteria_2 <https://doc.cgal.org/latest/Mesh_2/classCGAL_1_1Delaunay__mesh__criteria__2.html>`_
- `Delaunay_mesher_2 <https://doc.cgal.org/latest/Mesh_2/classCGAL_1_1Delaunay__mesher__2.html>`_
  
Attributes:
- **constraints** - Edges in the plane that encloses the mesh.
- **cdt** - `Constrained_triangulation_plus_2 <https://doc.cgal.org/latest/Triangulation_2/classCGAL_1_1Constrained__triangulation__plus__2.html>`_ object.
- **plane** -`CGAL Plan_3 <https://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Plane__3.html>`_ object.
- **edges** - Stores edges between subdomains with tags, even after removal of subdomain.

)doc";


/*static const char *__doc_Slice_get_points    = R"doc(Returns points in the 2D triangulation as a 2D numpy array.)doc";

static const char *__doc_Slice_get_facets_points = R"doc(Returns the points for facets in the 2D triangulation.)doc";

static const char *__doc_Slice_get_facet_tags = R"doc(Returns the points for facets in the 2D triangulation.)doc";

static const char *__doc_Slice_get_edge_tags = R"doc(Returns the points for facets in the 2D triangulation.)doc";

static const char *__doc_Slice_get_edges = R"doc(Returns the points for facets in the 2D triangulation.)doc";

static const char *__doc_Slice_get_edges_2 = R"doc(Returns the edges with the edge index. :Returns: Dictonary of edges and edge index. )doc";
*/ 



static const char *__doc_Slice_get_points =
R"doc( Returns points in the 2D triangulation as a numpy array with dim Nx3

:Returns: points in the mesh.

)doc";

static const char *__doc_Slice_get_cell_tags =
R"doc( Returns a numpy array of cell tags corresponding to the ouput of :func:`Slice.get_cells`. 


:Returns: numpy array of integers, represented as cell tags.

)doc";

static const char *__doc_Slice_get_cells =
R"doc( Returns all cells in the 2D triangulation, with a cell given by three vertex indices, as a numpy array.

Note : Vertex indices starts with index 1.

:Returns: numpy array of cells in the 2D triangulation.

)doc";


static const char *__doc_Slice_get_facets =
R"doc( Returns all facets in the 2D triangulation, with a facet given by two vertex indices, as numpy array, with the option of excluding unmarked facets.

Note : Vertex indices starts with index 1.

:param exclude_unmarked: if true exlcudes unmarked facets from output. 

:Returns: numpy array of facets in the mesh.

)doc";

static const char *__doc_Slice_get_facet_tags =
R"doc( Returns a numpy array of facet tags corresponding to the ouput of :func:`Slice.get_facets`. 

:param exclude_unmarked: if true exlcudes unmarked facets from output. 

:Returns: numpy array of facet tags in the mesh.


)doc";


static const char *__doc_Slice_add_polygon_domain = 
R"doc( Marks all cells within the polygon with the polygon tag.

:param polygon: a closed list of :class:`Point_2` objects.

:param polygon_tag: a integer for the polygon.

)doc";

static const char *__doc_Slice_Slice = R"doc(Create an empty SVMTK Slice object.)doc";

static const char *__doc_Slice_Slice_2 =
R"doc(Constructs a SVMTK Slice object with a specified plane. 

:param plane: Plane of the slice.

)doc";

static const char *__doc_Slice_Slice_3 =
R"doc(Constructs a SVMTK Slice object with a specified plane. 

:param point: :class:`Point_3` object on plane of the slice.
:param vector: :class:`Vector_3` plane normal.

)doc";

static const char *__doc_Slice_Slice_4 =
R"doc(Constructs a SVMTK Slice object with a specified plane_3. The plane is defined with the equation 
.. :math:`ax+by+cz+d = 0`

:param a: parameter in the plane equation.
:param b: parameter in the plane equation.
:param c: parameter in the plane equation.
:param d: parameter in the plane equation.

)doc";
    
static const char *__doc_Slice_remove_polygons = 
R"doc(Removes small closed loops of points, i.e. polygons that are lower than an area threshold

:param area_threshold: lower bound of the polygon area.

)doc";    

static const char *__doc_Slice_add_constraint =
R"doc(Adds polyline to constraints .

:param polyline: List of sequential :class:`Point_2` objects.

)doc";

static const char *__doc_Slice_add_constraints =
R"doc(Adds constraints from another Slice object.

:param slice: :class:`Slice` object.

)doc";

static const char *__doc_Slice_add_constraints_2 =
R"doc(Adds polylines to constraints.

:param polylines: Nested list of sequential :class:`Point_2` objects.

)doc";

static const char *__doc_Slice_add_surface_domains =
R"doc(Add tags to the facets in the 2D mesh based on overlapping surfaces and :class:`SubdomainMaps`.

Adds tags to faces dependent on position (inside/outside) according to closed triangulated surface in 3D.

:param surfaces: List of :class:`Surface` objects
:param map: :class:`SubdomainMaps` object. 

)doc";

static const char *__doc_Slice_add_surface_domains_2 =
R"doc(Add tags to the facets in the 2D mesh based on overlapping surfaces and DefaultMap.

:param surfaces: List of :class:`Surface` objects.

)doc";

static const char *__doc_Slice_clear_costraints = R"doc(Clears all the constraints added to the class object)doc";

static const char *__doc_Slice_connected_components =
R"doc(Returns the number of connected components.

:Returns: The number of connected components.

)doc";

static const char *__doc_Slice_create_mesh =
R"doc(Create 2D mesh given a set of constraints, i.e. specifiec edges.

:param mesh_resolution: The value divisor used with the minimum bounding radius to determine the maximum edge size. 

)doc";


static const char *__doc_Slice_create_mesh_2 =
R"doc(Create 2D mesh given a set of constraints, i.e. specifiec edges.

:param min_angle: Sets the lower-bound for the angles (in degrees) of the mesh facets.
:param edge_size: Sets the lower-bound of the edge size of the mesh.
)doc";


static const char *__doc_Slice_export_as_surface =
R"doc(Transforms the 2D mesh to 3D surface mesh stored in a SVMTK Surface
object. 

:Returns: :class:`Surface` object.

:Raises: *EmptyMeshError* if mesh variable is empty.

)doc";

static const char *__doc_Slice_get_constraints =
R"doc(Get the constraints added to the class object

:Returns: The constraints added to the class object.

)doc";



static const char *__doc_Slice_get_plane =
R"doc(Get the plane attribut.

:Returns: :class:`Plane_3` of the object.

)doc";

static const char *__doc_Slice_get_subdomains =
R"doc(Returns the subdomain tags in the mesh.

:Returns: Set of integers that represents the subdomain tags in the mesh.

:Raises: *EmptyMeshError* if mesh variable is empty.

)doc";

static const char *__doc_Slice_keep_largest_connected_component = R"doc(Calculates and keeps the largest connected component.)doc";

static const char *__doc_Slice_num_constraints =
R"doc(Returns the number of constraints.

:Returns: The number of constraints.

)doc";

static const char *__doc_Slice_num_cells =
R"doc(Returns the number of cells.

:Returns: The number of cells, i.e triangles.

)doc";

static const char *__doc_Slice_num_facets =
R"doc(Returns the number of facets.

:Returns: The number of facets, i.e edges.

)doc";

static const char *__doc_Slice_num_vertices =
R"doc(Returns the number of vertices.

:Returns: The number of vertices.

)doc";

static const char *__doc_Slice_num_subdomains =
R"doc(Returns the number of subdomains.

:Returns: The number of subdomains.

)doc";

static const char *__doc_Slice_remove_subdomain =
R"doc(Removes all faces in the mesh with a specified tag.

:param tag: removes facets with this integer tag.

)doc";

static const char *__doc_Slice_remove_subdomain_2 =
R"doc(Removes all faces in the mesh with specificed tags.

:param tags:  List of facet tags to be removed.

:Raises: *EmptyMeshError* if mesh variable is empty.

)doc";

static const char *__doc_Slice_save =
R"doc(Saves the 2D mesh to file. Valid formats : off, stl, vtu and mesh (with tags).

:param filename: filename of the 2D mesh is to be stored.

)doc";

static const char *__doc_Slice_set_plane =
R"doc(Set the plane attribute.

:param plane: :class:`Plane_3` object. 

)doc";

static const char *__doc_Slice_simplify =
R"doc(Simplify number of constraints to be a specific fraction.

:param threshold: Stops when ratio of remaining vertices is lower than the threshold. 

)doc";

static const char *__doc_Slice_slice_surfaces =
R"doc(Slices surfaces with the :class:`Slice` plane and adds the intersections as constraints in object.

:param surfaces: List of :class:`Surface` objects.

)doc";

static const char *__doc_Slice_slice_mesh =
R"doc(Slices a volumetric mesh with the :class:`Slice` plane and adds the intersections as constraints in object.

:param mesh: :class:`Domain` object.

)doc";


static const char *__doc_Slice_bifurcation_split =
R"doc( Divides constraints with a bifurcation point into separate constraints 

If a point is shared with multiple grahps, then split all edge to point into different constraints.


)doc";


static const char *__doc_SubdomainMap =
R"doc(This class is used for setting the subdomain tag during the construction of the mesh. 


A vertex position relative to the surface can be written as a bitstring, i.e. inside or outside of the surface. The evaluation of the binary string that can be used to map a point to have a specific tag.

In the construction of the mesh, cells will be tagged based on the location of the cell relative to the surfaces ovelaps.


Attributes: 
- **num_surfaces** - stores the number surfaces.   
- **subdmap** - subdomain map, bitset integer dictonary.
- **patches** - interface map stored as dictonary with integer tuple as key and i integermap.


)doc";



static const char *__doc_SubdomainMap_SubdomainMap =
R"doc( Creates SVMKT SubdomainMap object.

:param num_surfaces: The number of surfaces. Default value = 0.

)doc";

static const char *__doc_SubdomainMap_add =
R"doc(Transforms a string of 0 and 1 to a binarystring.

:param string: bitstring, i.e. contains only 0 and 1.
:param tag: The tag that the specified subdomain will have in the stored mesh.

)doc";

static const char *__doc_SubdomainMap_add_interface =
R"doc(Adds a tag value for surfaces patches between subdomains defined by a pair of integer the class attribute patches.

:param interface: Tuple of two integers that defines the surface interface between two subdomains.
:param tag: The tag that the specified surface interface will have in the stored mesh.

)doc";

static const char *__doc_SubdomainMap_erase =
R"doc(Erase binary string from SubdomainMap.

:param string: Removal of bitstring from SubdomainMap.

)doc";

static const char *__doc_SubdomainMap_get_interfaces =
R"doc(Get stored interface patches and tags.

)doc";

static const char *__doc_SubdomainMap_get_map =
R"doc(Returns the map.

:Returns: List of integer that represents the tags added to the class object.

)doc";

static const char *__doc_SubdomainMap_get_tags =
R"doc(Returns all tags that is added to the class object.

:Returns: List of integer that represents the tags added to the class object.

)doc";


static const char *__doc_SubdomainMap_make_interfaces =
R"doc(Returns the content of class attribute patches between subdomains with the corresponding tag value. If patches is empty, gives each interfaces an unique tag based on presence in mesh and the highest cell tag.

:param interfaces: Tuple of two integers that defines the suface interface between two subdomains.

:Returns: Dictonary with tuple of two integers as key and a integer tag value.

)doc";


static const char *__doc_SubdomainMap_print =
R"doc(Prints subdomain binarystring and patches map.)doc";

static const char *__doc_SubdomainMap_set_num_surfaces =
R"doc(Sets the number of surfaces

The number of surfaces is used to automatically fill binary combinations when adding a sting with asterix. 

:param number_of_surfaces: An integer equal the number of surfaces.
    
)doc";

static const char *__doc_Surface =
R"doc(The SVMTK Surface class is used to create, store and manipulate triangulated surfaces in 3D. 

CGAL is implemented with different kernels, that have different properties, and this class is implemnented with the `Exact predicates inexact constructions kernel <https://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Exact__predicates__inexact__constructions__kernel.html>`_.

SVMTK Surface class is used to handle operations related to triangualted surfaces in 3D CGAL `Surface_mesh <https://doc.cgal.org/latest/Surface_mesh/index.html>`_. It should be noted that most operations can also be used with CGAL `Polyhedron_3 <https://doc.cgal.org/latest/Polyhedron/index.html>`_.

Notable CGAL declarations: 
- `Surface_mesh_default_triangulation_3 <https://doc.cgal.org/latest/Surface_mesher/classCGAL_1_1Surface__mesh__default__triangulation__3.html>`_.
- `Side_of_triangle_mesh <https://doc.cgal.org/latest/Polygon_mesh_processing/classCGAL_1_1Side__of__triangle__mesh.html>`_.
 
Attributes:
- **mesh** - `Surface_mesh <https://doc.cgal.org/latest/Surface_mesh/index.html>`_ object.

)doc";

static const char *__doc_Surface_Surface = R"doc(Constructs an empty SVMTK Surface object.)doc";

static const char *__doc_Surface_Surface_2 =
R"doc(Constructs a SVMTK Surface object with surface from file. Current supported fileformats: off and stl.

:parameters filename: The filename of the surface to be loaded into the :class:`Surface` object.

)doc";

static const char *__doc_Surface_Surface_3 = R"doc(Constructs a copy of a SVMTK Surface object.)doc";

static const char *__doc_Surface_adjust_boundary =
R"doc(Moves the surfaces vertices in the normal vertex direction multiplied with a specified value.

:param c: Multiplier of the unit normal vertex vector.

)doc";

static const char *__doc_Surface_area =
R"doc(Computes the area of the surface.

:Returns: The area of the surface.

)doc";

static const char *__doc_Surface_average_edge_length =
R"doc(Computes the average edge length in the stored mesh object.

:Returns: The average edge length.

)doc";

static const char *__doc_Surface_centeroid =
R"doc(Computes the centeroid point of the surface.

:Returns: :class:`Point_3` object of the surface centeroid.


)doc";

static const char *__doc_Surface_clear = R"doc(Clears mesh.)doc";

static const char *__doc_Surface_copy = R"doc(Creates a deepcopy of the :class:`Surface` object.)doc";


static const char *__doc_Surface_clip =
R"doc(Clips the surface mesh for a given :class:`Plane_3` object.

This will construct a circle that spans the surface mesh in a given plane and the circle center is the projection of the surface centeroid onto the plane. Then, clipping ths surface mesh with the constructed surface. This method will cause the faces on the clipped boundary to have similar size as the rest of the mesh. The clipping is done using the CGAL function `clip <https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__corefinement__grp.html#ga30082762ba2d947cba304e2884d96a99>`_.

:param a: parameter in plane equation.
:param b: parameter in plane equation.
:param c: parameter in plane equation.
:param d: parameter in plane equation.
:param preserve_manifold: True to preserve manifold.

:Returns: True if successful.

)doc";

static const char *__doc_Surface_clip_2 =
R"doc(Clips the surface mesh for a given :class:`Plane_3` object.

:param point: :class:`Point_3` on the plane.
:param vector: :class:`Vector_3` object normal to the plane.
:param preserve_manifold: - true to preserve manifold.

:Returns: True if successful.

)doc";

static const char *__doc_Surface_clip_3 =
R"doc(Clips the surface mesh for a given :class:`Plane_3` object.

:param plane: :class:`Plane_3` object.
:param preserve_manifold:  True will preserve manifold.

:Returns: True if successful.


)doc";

static const char *__doc_Surface_clip_4 =
R"doc(Clips a surface mesh given another SVMTK Surface object.

:param other: :class:`Surface` object.
:param preserve_manifold: True to preserve manifold.
:param invert: True will invert clip plane.

:Returns: True if successful.

)doc";

static const char *__doc_Surface_clip_5 =
R"doc(Clips the surface mesh with circle in a specified plane.

Creates a circle in a specified plane, and uses this circle to clip the surface mesh. The clips made using a surface mesh creates better triangulation on the clipped boundary, and can be used to perform a localized clip.

:param point: :class:`Point_3` on the plane.
:param vector: :class:`Vector_3` object normal to the plane.
:param radius: Radius of the circle.
:param invert: True will invert clip plane.
:param preserve_manifold: True will preserve manifold

:Returns: True if successful.

)doc";


static const char *__doc_Surface_get_perpendicular_cut=
R"doc(Constructs a circular surface in a given plane, which is intersects the perpendicular to the :func:`mean_curvature_flow.

Computes the nearest point in the mean curvature flow to the point query, and creates a circular surface in the plane that is perpendicular to mean curvature flow.

See:
`mean_curvature_flow <https://doc.cgal.org/latest/Surface_mesh_skeletonization/index.html>`_.

:param point: :class:`Point_3` query to located the nearest point of the mean curvature flow.
:param radius: Radius of the circle.

:Returns: :class:`Surface` 

)doc";

static const char *__doc_Surface_get_perpendicular_cut_2=
R"doc(Constructs a circular surface in a given plane, which is intersects the perpendicular to the :func:`mean_curvature_flow.

Computes where the plane and mean curvature flow intersectes, and creates a circular surface in the plane that is perpendicular to mean curvature flow.

See:
`mean_curvature_flow <https://doc.cgal.org/latest/Surface_mesh_skeletonization/index.html>`_.

:param plane: :class:`Plane_3`. 
:param radius: Radius of the circle.

:Returns: :class:`Surface` that can be used to clip the surface.

)doc";


static const char *__doc_Surface_collapse_edges =
R"doc(Combines smaller edges together, so that all edges are larger than the input parameter

See also:
`edge_collapse <https://doc.cgal.org/latest/Surface_mesh_simplification/index.html>`_.


:param target_edge_length: The edge length that is targeted in the mesh by combining smaller edges. 


:Returns: The number of collapsed edges.

)doc";

static const char *__doc_Surface_collapse_edges_2 =
R"doc(Collapses smaller edges together if these edges can be represented as a longer edge.

:Returns: The number of collapsed edges.

)doc";


static const char *__doc_Surface_remove_degenerate_faces =
R"doc( Finds and removes degenerate faces in the surface and subsequent hole filling of the surface.

)doc";

static const char *__doc_Surface_make_circle_in_plane_2 =
R"doc( Finds and removes degenerate faces in the surface.

:param point: the center point of the circle. 
:param vec  : the normal of the plane. 
:param radius : the radius of the circle. 
:param edge_length : the target edge length for the surface construction.  


)doc";







static const char *__doc_Surface_convex_hull =
R"doc(Computes the convex hull of the surface points.

See also:
`convex_hull_3 <https://doc.cgal.org/latest/Convex_hull_3/index.html>`_.

:Returns: :class:`Surface` object.

)doc";
    
static const char *__doc_Surface_cylindrical_connection =
R"doc(Constructs a cylindrical connection bridge between the shortests line between two points in different surface.

The user can use the union operation to combine it with the main mesh.

:param other: :class:`Surface` object.
:param radius: radius of the cylinder surface.
:param edge_length: The edge length that is targeted in the mesh. 

:Returns: :class:`Surface` object.

)doc";

static const char *__doc_Surface_cylindrical_extension =
R"doc(Constructs a cylindrical extension, which is a cylinder surface mesh combined with a sphere on the end closest to the mesh.The user can use the union operation to combine it with the main mesh. The cylinder surface mesh is determined by centeriod of vertex points that are closest to a point outisde the mesh, radius, length and the option of creating a cylinder normal to the surface mesh. Works best on convex surfaces.

:param point: :class:`Point_3` object outside the surface mesh, the closest surface vertex to point is found and used as an initial point.  
:param radius: Radius of the cylinder surface.
:param length: Length of the cylinder surface from the surface mesh.
:param edge_length: The edge length that is targeted in the mesh. 
:param normal: True will set the cylinder to be normal to the surface mesh, otherwise it will be in the direction of the initial point.

:Returns: :class:`Surface` object.

:Raises: *InvalidArgumentError* if point is inside surface.

)doc";

static const char *__doc_Surface_cylindrical_extension_2 =
R"doc(

:param x: x-coordinate of point.
:param y: y-coordinate of point.
:param z: z-coordinate of point.
:param radius: Radius of the cylinder.
:param length: Length of the cylinder.
:param edge_length: The edge length that is targeted in the mesh. 
:param normal: True will set the cylinder to be normal to the surface mesh, otherwise it will be in the direction of the initial point.

:Returns: :class:`Surface` object.

)doc";

static const char *__doc_Surface_distance_to_point =
R"doc(Computes the distance between a point and the closest vertex in the triangulated surface.

:param point: :class:`Point_3` object.

:Returns: The distance from point and the closest vertex on the surface.

)doc";

static const char *__doc_Surface_embed =
R"doc(Moves incrementally vertices in a negative normal direction so that they are inside specified surface.

:param other: :class:`Surface` object.
:param adjustment:  Multiplier of the smallest connected edge of a vertex, which equals the displacement of the vertices for each iteration. 
:param smoothing: Laplacian smoothing factor for each increment, see :func:`laplacian_smoothing`.
:param max_iter: The maximum number of iterations for a while loop.

:Returns: True if completed and the number vertices still outside surface.

)doc";

static const char *__doc_Surface_enclose =
R"doc(Moves incrementally vertices in a positiv normal direction so that they are outside specified surface.

:Precondition: The surfaces must not intersect.

:param other: :class:`Surface` object.
:param adjustment:  Multiplier of the smallest connected edge of a vertex, which equals the displacement of the vertices for each iteration. 
:param smoothing: Laplacian smoothing factor for each increment, see :func:`laplacian_smoothing`.
:param max_iter: The maximum number of iterations for a while loop.

:Returns: True if completed and the number vertices still inside surface.


)doc";

static const char *__doc_Surface_expose =
R"doc(Moves incrementally vertices in a negative normal direction so that they are outside specified surface.


:param other: :class:`Surface` object.
:param adjustment: Multiplier of the smallest connected edge of a vertex, which equals the displacement of the vertices for each iteration. 
:param smoothing: Laplacian smoothing factor for each increment, see :func:`laplacian_smoothing`.
:param max_iter: The maximum number of iterations for a while loop.

:Returns: True if completed and the number vertices still outside surface. 


)doc";

static const char *__doc_Surface_set_proximity_ratio =
R"doc(Sets the proximity ratio. 

Note: If the distance between two vertices exceeds the proximity ratio times edge length,
      then the vertices are considered close. Value should be between 1 and 2. 

:param proximity_ratio: used to determine close vertice:  

)doc";

static const char *__doc_Surface_set_smoothing_reduction =
R"doc( Sets the smoothing_reduction. 

Note: When the displacment is lower than the lower displacment bound, 
      then the smoothing reduction is used to decrease the smoothing.  
      
:param smooth_reduction: used to decrease smoothing if incremental displacement is lower than the bound. 


)doc";
static const char *__doc_Surface_set_displacment_ratio =
R"doc(Sets the displacement ratio. 

Note: The displacment ratio times average egde length determines the lower displacement bound 

:param displacment_ratio: 

)doc";




static const char *__doc_Surface_fill_holes =
R"doc(Finds and fills holes in surface mesh. 

Uses CGAL function `triangulate_refine_and_fair_hole <https://doc.cgal.org/latest/Polygon_mesh_processing/group__hole__filling__grp.html>`_ and follows `example <https://doc.cgal.org/latest/Polygon_mesh_processing/Polygon_mesh_processing_2hole_filling_example_SM_8cpp-example.html>`_.

:Returns: The number of holes filled.

)doc";

static const char *__doc_Surface_get_closest_points =
R"doc(Finds a specified number mesh points that are closest to a point not on the mesh.

:param source: :class:`Point_3` object not on the surface mesh. 
:param num: The number of points to return. 

:Returns: List of :class:`Point_3`.

)doc";

static const char *__doc_Surface_get_points =
R"doc(Returns the points of the surfaces mesh.

:Returns: List of :class:`Point_3`.

)doc";

static const char *__doc_Surface_get_shortest_surface_path =
R"doc(Computes and returns the shortest surface path between two points. Two points are projected on to the surface mesh, and the shortest surface path is found.

:param source: :class:`Point_3` object.
:param target: :class:`Point_3` object.

:Returns: List of sequential :class:`Point_3` objects.

)doc";
    
static const char *__doc_Surface_get_shortest_surface_path_2 =
R"doc(Computes and returns the shortest surface path between two points.

:param x0: x-coordinate of the source point. 
:param y0: y-coordinate of the source point. 
:param z0: z-coordinate of the source point. 
:param x1: x-coordinate of the target point. 
:param y1: y-coordinate of the target point. 
:param z1: z-coordinate of the target point. 

:Returns: List of sequential :class:`Point_3` objects that is the shortest surface path.

)doc";

static const char *__doc_Surface_repair_self_intersections =
R"doc(Removes self intersection from surface. This function uses experimental CGAL algorithms to remove self-intersections, see `remove_self_intersections <https://github.com/CGAL/cgal/blob/master/Polygon_mesh_processing/include/CGAL/Polygon_mesh_processing/repair_self_intersections.h>`_.

:param volume_threshold: Ratio value of the volume such that only connected components whose volume is larger than this value are kept (only applies to closed connected components)
:param cap_thresold: Anlge between 160 180 degrees CGAL explanation->the cosine of a minimum angle such that if `f` has an angle greater than this bound, it is a cap. The threshold is in range `[-1 0]` and corresponds to an angle between `90` and `180` degrees.
:param needle_threshold: Long edge divided by the short edge of a cell.
:param collapse_threshold: Ratio of the average edge length.

:Returns: True if complete and number of remaining self-intersections.

)doc";
    
static const char *__doc_Surface_get_slice =
R"doc(Slices a SVMTK :class:`Surface` object according to a given plane, and returns a SVMTK Slice object.

:param plane: :class:`Plane_3` object used to slice the :class:`Surface` object. 

:Returns: :class:`Slice` object.

)doc";

static const char *__doc_Surface_get_slice_2 =
R"doc(Slices a surface mesh based on a plane, that is defined by the plane equation :math:`ax+by+cz+d = 0`. Better result if the plane normal is close to unity.

:param a: plane equation parameter. 
:param b: plane equation parameter.
:param c: plane equation parameter.
:param d: plane equation parameter.

:Returns: :class:`Slice` object.

)doc";

static const char *__doc_Surface_implicit_surface =
R"doc(Creates a surface mesh based on an implicit function

See also:
`Surface_mesher <https://doc.cgal.org/latest/Surface_mesher/index.html>`_.


:param implitict_function: Python function defined as f(x,y,z)=0 and the interior defined as f(x,y,z) < 0.
:param bounding_sphere_radius: Radius of a sphere that encloses the mesh construction.
:param angular_bound: Bounds for the minimum facet angle in degrees.
:param radius_bound: Bound for the minimum for the radius of the surface Delaunay balls and the center-center distances respectively.
:param distance_bound: Bound for the minimum center-center distances respectively.

)doc";

static const char *__doc_Surface_is_point_inside =
R"doc(Checks if a point is inside surface mesh.

Query if a point is inside the surface mesh.

:param point: :class:`Point_3` object.

:Returns: True if point is inside surface otherwise false.

)doc";

static const char *__doc_Surface_isotropic_remeshing =
R"doc(Isotropic remeshing of surface mesh. Remeshing of the surface mesh so that all edges have the same length.

Uses CGAL `split_long_edges <https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__meshing__grp.html>`_. and `isotropic_remeshing <https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__meshing__grp.html>`_. Split_long_edges avoids the pitfall described `here <https://doc.cgal.org/5.0.3/Polygon_mesh_processing/index.html>`_.

:param target_edge_length: The edge length that is targeted in the remeshed patch. If 0 is passed then only the edge-flip, tangential relaxation, and projection steps will be done.
:param nb_iter: the number of iterations for the sequence of atomic operations performed. 
:param protect_border: If true, constraint edges cannot be modified at all during the remeshing process.
    
)doc";

static const char *__doc_Surface_keep_largest_connected_component =
R"doc(Keeps the largest connected component of the stored mesh.

See also:
`keep_largest_connected_components <https://doc.cgal.org/latest/Polygon_mesh_processing/group__keep__connected__components__grp.html>`_.

:Returns: The number of connected components.

)doc";

static const char *__doc_Surface_connected_components =
R"doc(Creates :class:`Surface` object for each connected component.


:Returns: List of :class:`Surface` object for each connected component 

)doc";


static const char *__doc_Surface_make_circle_in_plane =
R"doc(Constructs a circle in the a given surface plane.


:param px: x-coordinate of circle center in plane.  
:param py: y-coordinate of circle center in plane. 
:param pz: z-coordinate of circle center in plane.  
:param vx: x-coordinate of the plane normal vector.
:param vy: y-coordinate of the plane normal vector.
:param vz: z-coordinate of the plane normal vector.
:param radius: radius of circle.
:param edge_length: The upper edge length of the constructed triangulated surface. 

)doc";

static const char *__doc_Surface_make_cone =
R"doc(Creates a surface mesh structure with vertices and facets connecting vertices for a cone. The function also handles the special cases of sharp cone and cylinder.

:Parameters: |
:param x0: x-coordinate of the first cone center.
:param y0: y-coordinate of the first cone center.
:param z0: z-coordinate of the first cone center.
:param x1: x-coordinate of the second cone center.
:param y1: y-coordinate of the second cone center.
:param z1: z-coordinate of the second cone center.
:param r0: radius corresponding to the first cone center.
:param r1: radius corresponding to the second cone center.
:param edge_length: The edge length that is targeted in the mesh. 

)doc";

static const char *__doc_Surface_make_cone_2 =
R"doc(Constructs a cone surface.

:param p0: :class:`Point_3` object of the bottom center of the cone.
:param p1: :class:`Point_3` object of the top center of the cone.
:param r0: the bottom radius of the cone.
:param r1: the top radius of the cone.
:param edge_length: The edge length that is targeted in the mesh. 

)doc";

static const char *__doc_Surface_make_cube =
R"doc(Creates a surface mesh structure with vertices and facets connecting vertices for a cube.

:param x0: Sets the x-coordinate of the first cube corner point. 
:param y0: Sets the y-coordinate of the first cube corner point. 
:param z0: Sets the z-coordinate of the first cube corner point. 
:param x1: Sets the x-coordinate of the second cube corner, opposite of the first corner. 
:param y1: Sets the y-coordinate of the second cube corner, opposite of the first corner.  
:param z1: Sets the z-coordinate of the second cube corner, opposite of the first corner.  
:param edge_length: The edge length that is targeted in the mesh. 

)doc";

static const char *__doc_Surface_make_cube_2 =
R"doc(Constructs a cube surface.

:Parameters: |
:param p0: :class:`Point_3` object that sets the first cube corner
:param p1: :class:`Point_3` object that sets the second cube corner, opposite of the first corner.
:param edge_length: The edge length that is targeted in the mesh. 

)doc";

static const char *__doc_Surface_make_cylinder =
R"doc(Creates a surface mesh structure with vertices and facets connecting vertices for a cylinder.

:Parameters: |
:param x0: Sets the x-coordinate of the bottom cylinder center.
:param y0: Sets the y-coordinate of the bottom cylinder center.
:param z0: Sets the z-coordinate of the bottom cylinder center.
:param x1: Sets the x-coordinate of the top cylinder center.
:param y1: Sets the y-coordinate of the top cylinder center.
:param z1: Sets the z-coordinate of the top cylinder center.
:param r0: Sets the radius of the cylinder.
:param edge_length: The edge length that is targeted in the mesh. 

)doc";

static const char *__doc_Surface_make_cylinder_2 =
R"doc(Creates a surface mesh structure with vertices and facets connecting vertices for a cylinder.

:param p0: :class:`Point_3` object of the bottom center of the cylinder.
:param p1: :class:`Point_3` object of the top center of the cylinder.
:param r0: Sets the radius of the cylinder.
:param edge_length: the upper edge length of the constructed triangulated surface.


)doc";

static const char *__doc_Surface_make_sphere =
R"doc(Creates a sphere surface mesh. Uses an implicit function to construct a structure of vertices and facets connecting vertices in the shape of a sphere.

:param x0: Sets x-coordinate of the sphere center.
:param y0: Sets y-coordinate of the sphere center.
:param z0: Sets z-coordinate of the sphere center.
:param r0: Sets the radius of the sphere.
:param mesh_resolution: ratio between the sphere radius divided by the maximum edge length of the resulting surface mesh.

)doc";

static const char *__doc_Surface_make_sphere_2 =
R"doc(Creates a sphere surface mesh.

:param point: :class:`Point_3` object of the sphere center. 
:param radius: Sets the radius of the sphere.
:param edge_length: The edge length that is targeted in the mesh.

)doc";

static const char *__doc_Surface_mean_curvature_flow =
R"doc(Computes the centerline of the surface.

See also:
`extract_mean_curvature_flow_skeleton <https://doc.cgal.org/latest/Surface_mesh_skeletonization/index.html>`_.

:Returns: List of sequential :class:`Point_3` objects.
    
)doc";

static const char *__doc_Surface_num_edges =
R"doc(Returns the number of edges in the surface.

:Return0s: The number of edges in the surface.

)doc";

static const char *__doc_Surface_num_faces =
R"doc(Returns the number of faces in the surface.

:Returns:  The number of faces in the surface.

)doc";

static const char *__doc_Surface_num_self_intersections =
R"doc(Returns number of self intersection in the surface.

:Returns: The number of self intersection in the surface.

)doc";

static const char *__doc_Surface_num_vertices =
R"doc(Returns the number of vertices in the surface.

:Returns: The number of vertices in the surface.

)doc";

static const char *__doc_Surface_reconstruct =
R"doc(Reconstruct surface mesh.

Reconstruct a surface based on a CGAL surface mesh object with points using CGAL `poisson_reconstruction <https://doc.cgal.org/latest/Poisson_surface_reconstruction_3/index.html>`_. 

:param bounding_sphere_radius: The radius of a sphere bounding the meshing operations
:param angular_bound: Bounds for the minimum facet angle in degrees.
:param radius_bound: Bound for the minimum for the radius of the surface Delaunay balls and the center-center distances respectively. 
:param distance_bound: Bound for the minimum center-center distances respectively.


)doc";

static const char *__doc_Surface_reconstruct_2 =
R"doc(Reconstruct surface mesh.

Reconstruct a surface based on a CGAL surface mesh object with points using CGAL `poisson_reconstruction <https://doc.cgal.org/latest/Poisson_surface_reconstruction_3/index.html>`_. 

:param filenam: Filename of the data file containing points with normals.
:param bounding_sphere_radius: The radius of a sphere bounding the meshing operations
:param angular_bound: Bounds for the minimum facet angle in degrees.
:param radius_bound: Bound for the minimum for the radius of the surface Delaunay balls and the center-center distances respectively. 
:param distance_bound: Bound for the minimum center-center distances respectively.


)doc";


static const char *__doc_Surface_remove_small_components =
R"doc(Removes connected components whose area or volume is under a certain threshold value.

:Precondition: The surface mesh must bound a volume.


:param volume_threshold: Ratio value of the volume such that only connected components whose volume is larger are kept (only applies to closed connected components).


Raises: *PreconditionError* if surface does not enclose volume.


See also:
`remove_connected_components_of_negligible_size <https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__repairing__grp.html#gac544fcaba1d59d330a3a1536caff392a>`_.

)doc";



static const char *__doc_Surface_save =
R"doc(Saves the surface mesh to file. Valid file formats: off and stl.

:param outpath: A string path to save file. 

)doc";

static const char *__doc_Surface_separate =
R"doc(Moves incrementally vertices that are close to another surface in opposite direction. 

:Precondition: The other surface must either be enclosed or embedded in the surface.

:Parameters:
:param other: :class:`Surface` object.
:param adjustment:  Multiplier of the smallest connected edge of a vertex, which equals the displacement of the vertices for each iteration. 
:param smoothing: Laplacian smoothing factor for each increment, see :func:`laplacian_smoothing`.
:param max_iter: The maximum number of iterations. 

:Returns: True if completed and number of remaining vertices to move.

:Raises: *InvalidArgumentError* if surfaces intersect.


)doc";

static const char *__doc_Surface_separate_close_vertices =
R"doc(Incrementally separates close non-connected surface mesh vertices.

Separates close non-connected surface mesh vertices so that the closest surface mesh vertex is one that is connected, i.e. shares a edge.

:param adjustment:  Multiplier of the smallest connected edge of a vertex, which equals the displacement of the vertices for each iteration. 
:param max_iter: the maximum number of iterations.

:Returns: True if completed and number of remaining vertices to move.

)doc";

static const char *__doc_Surface_separate_narrow_gaps =
R"doc(Separates unconnected surface mesh vertices that are close in a positive normal direction by moving vertices in a negative surface normal direction, i.e. contraction.

:param adjustment: Multiplier of the smallest connected edge of a vertex, which equals the displacement of the vertices for each iteration. 

:Returns: True if completed and number of remaining vertices to move.

)doc";

static const char *__doc_Surface_smooth_laplacian =
R"doc(Smooths all vertices of surface mesh with Laplacian smoothing.

This is done by taking the sum of the vector edges for each vertex, and multiplied with the constant float input.

:param c: multipler of the displacement vector that gives the new vertex coordinates.
:param nb_iter: The number of iterations of Laplacian smoothing

)doc";

static const char *__doc_Surface_smooth_shape =
R"doc(Smooth the triangulated surface.

See also:
`smooth_shape <https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__meshing__grp.html>`_.

:param time: step that corresponds to the speed by which the surface is smoothed. A larger time step results in faster convergence but details may be distorted to have a larger extent compared to more iterations with a smaller step. Typical values scale in the interval (1e-6, 1].
:param nb_iter: the number of iteations. 

)doc";

static const char *__doc_Surface_smooth_taubin =
R"doc(Taubin smothing of surface mesh.

Taubin smoothing of the surface vertices. This corresponds to a Laplacian smoothing with value \lambda, followed by a Laplacian smoothing with value \mu Given the requirement: :math:`\lambda < -\mu`. The Laplacian smoothing parameters are set, but the user may construct with their own parameters with smooth_laplacian function.

:param nb_iter: The number of iterations of Taubin smoothing. 

)doc";

static const char *__doc_Surface_span =
R"doc(Finds the span in a given Cartesian direction, i.e. (min, max) vertex
position.

:param direction: Cartesian direction (0=x, 1=y, 2=z).

:Returns: The minimum coordinate and max coordinate value in a given direction.

)doc";

static const char *__doc_Surface_split_edges =
R"doc(splits edges to a target edge length CGAL function for splitting long
edges.

See also:
`split_long_edges <https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__meshing__grp.html>`_.

:param target_edge_length: The target edge length that the longer edges splits.

)doc";

static const char *__doc_Surface_surface_difference =
R"doc(Computes the difference between two triangulated surface mesh.

See also:
`corefine_and_compute_difference <https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__corefinement__grp.html>`_. 
    
:Precondition: Both surfaces bounds a volume and does not have self-intersections where the surfaces intersect.

:param other: :class:`Surface` object.

:Returns: True if intersection computation is successful.

)doc";

static const char *__doc_Surface_surface_intersection =
R"doc(Computes the intersection between two triangulated surface mesh.
  
The functions have precondition that both surfaces bounds a volume and that both surfaces does not have self-intersections where the surfaces intersect.

See also:
`corefine_and_compute_intersection <https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__corefinement__grp.html>`_. 
 
:Precondition: Both surfaces bounds a volume and does not have self-intersections where the surfaces intersect.

:param other: :class:`Surface` object.

:Returns: True if intersection computation is successful.

)doc";

static const char *__doc_Surface_surface_union =
R"doc(Computes the union between two triangulated surface mesh.

See also: `corefine_and_compute_union <https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__corefinement__grp.html>`_.
    
:Precondition: Both surfaces bounds a volume and does not have self-intersections where the surfaces intersect.

:param other: :class:`Surface` object.

:Returns: True if intersection computation is successful.

)doc";

static const char *__doc_Surface_triangulate_faces =
R"doc(Triangulates faces of the surface mesh. Uses CGAL function `triangulate_faces <https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__meshing__grp.html>`_.

:Returns: True if complete.

)doc";

static const char *__doc_Surface_volume =
R"doc(Computes the volume enclosed by the surface.

:Returns: The volume of the surface.
    
:Raises: *PreconditionError* if surface does not enclose volume.

)doc";


static const char *__doc_Vector3 =
R"doc(Wrapper for `CGAL Vector_3 class <https://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Vector__3.html>`_.

)doc";

static const char *__doc_Vector3_Vector3 =
R"doc(Creates a :class`Vector_3` object.

:param x: Sets x direction of the object.
:param y: Sets y direction of the object.
:param z: Sets z direction of the object. 

)doc";

static const char *__doc_Vector3_Vector3_2 =
R"doc(Creates a :class`Vector_3` object as difference between two points.

:param source: :class:`Point_3` object. 
:param target: :class:`Point_3` object.

)doc";


static const char *__doc_Vector2 =
R"doc(Wrapper for `CGAL Vector_2 class <https://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Vector__2.html>`_.

)doc";

static const char *__doc_Vector2_Vector2 =
R"doc(Creates a :class`Vector_2` object.

:param x: Sets x direction of the object.
:param y: Sets y direction of the object.


)doc";

static const char *__doc_Vector2_Vector2_2 =
R"doc(Creates a :class`Vector_2` object as difference between two points.

:param source: :class:`Point_2` object. 
:param target: :class:`Point_2` object.

)doc";


static const char *__doc_convex_hull =
R"doc(Creates the convex hull of the surface points and return :class:`Surface`. 

See `CGAL convex hull <https://doc.cgal.org/latest/Convex_hull_3/index.html>`_.

:Return: :class:`Surface_3`
)doc";

static const char *__doc_smooth_polylines =
R"doc( Smooth polyline in 2D, with the endpoints remain fixed.

Auxillary function to smooth a list of :class:`Point_2`. 

:param polyline: a list of :class:`Point_2`.
:param beta: the laplacian smoothin factor. 

:Return: smoothed list of :class:`Point_2`.
)doc";

static const char *__doc_smooth_polylines_2 =
R"doc(Smooth polyline in 3D, with the endpoints remain fixed.

Auxillary function to smooth a list of :class:`Point_3`. 

:param polyline: a list of :class:`Point_3`.
:param beta: the laplacian smoothin factor. 

:Return: smoothed list of :class:`Point_3`.

)doc";

static const char *__doc_embed =
R"doc(Embeds the primary surface inside the secondary, and subsequently separates the surfaces. 

The function will embed the primary surface into the secondary by incrementally shrink the surface until all 
surface points are inside the secondary surface. Then, the surfaces are separate until the closest point for all 
points in the primary surfaces are connected with an edge, i.e. not a point in the secondary surface. 
This function make use of the function :func:`Surface::embed` and :func:`Surface::separate`.


:param surf1: the primary :class:`Surface`. 
:param surf2: the secondary :class:`Surface`. 
:param adjustment: edge factor of the displacment for each iteration.  
:param smoothing: smoothing factor of the displacment for each iteration.  
:param max_iter: maximum number of iterations.  

:Returns: True if successful.

)doc";

static const char *__doc_expose =
R"doc(Disjoins two surfaces, and subsequently separates the surfaces. 

The function will disjoin two surfaces by incrementally shrink the primary surfaces until all 
surface points are outside the secondary surface. Then, the surfaces are separate until the closest 
point in the primary surfaces is connected with an edge, i.e. not a point in the secondary surface. 
This function make use of the function :func:`Surface::expose` and :func:`Surface::separate`


:param surf1: the primary :class:`Surface`. 
:param surf2: the secondary :class:`Surface`. 
:param adjustment: edge factor of the displacment for each iteration.  
:param smoothing: smoothing factor of the displacment for each iteration.  
:param max_iter: maximum number of iterations.  

:Returns: True if successful.

)doc";


static const char *__doc_enclose =
R"doc(Enclose the primary surface inside the secondary, and subsequently separates the surfaces. 

The function will cause the primary surface to enclose the secondary by incrementally expanding the surface until all 
surface points are outside the secondary surface. Then, the surfaces are separate until the closest point for all 
points in the primary surfaces are connected with an edge, i.e. not a point in the secondary surface.  
This function make use of the function :func:`Surface::embed` and :func:`Surface::separate`.


:param surf1: the primary :class:`Surface`. 
:param surf2: the secondary :class:`Surface`. 
:param adjustment: edge factor of the displacment for each iteration.  
:param smoothing: smoothing factor of the displacment for each iteration.  
:param max_iter: maximum number of iterations.  

:Returns: True if successful.

)doc";

static const char *__doc_inside_polygon =
R"doc( Checks if a point in 2D is inside a polygon.

Checks if a :class:`Point_2` is inside a closed list of :class:`Point_2`

:param polygon: a closed loop of :class:`Point_2`
:param query: :class:`Point_2` to be checked 

:Return: true if point is inside polygon otherwise false.

)doc";




static const char *__doc_separate_close_surfaces =
R"doc(Separates two close surfaces outside a third surface. 

Separates the points of two surfaces that are outside a third surface iteratively by moving the vertices in the negative normal direction that is determined by the multiplication of a negative value and the shortest edge length corresponding to each vertex.This continues until the closest vertex is adjacent.

:param surf1: :class:`Surface` object.
:param surf2: :class:`Surface` object.
:param other: :class:`Surface` object, the algorithm is not applied to vertices inside this surface.
:param edge_movement: Multiplier of the vector, vertex to the other surface, that gives the displacement of the vertices for each iteration. The length of the vector is less than the shortest connected edge of the vertex.
:param smoothing: Laplacian smoothing factor for each increment, see :func:`laplacian_smoothing`.
:param max_iter: The maximum number of iteration.

:Returns: True if complete.

)doc";

static const char *__doc_separate_close_surfaces_2 =
R"doc(Separates two close surfaces.

Separates two surfaces iteratively by moving the vertices in the negative normal direction that is determined by the multiplication of a negative value and the shortest edge length corresponding to each vertex. This continues until the closest vertex is adjacent.

:param surf1: :class:`Surface` object.
:param surf2: :class:`Surface` object.
:param edge_movement: Multiplier of the vector, vertex to the other surface, which equals the displacement of the vertices for each iteration. The length of the vector is less than the shortest connected edge of the vertex.
:param smoothing: Laplacian smoothing factor for each increment, see :func:`laplacian_smoothing`.
:param max_iter: The maximum number of iterations. 

:Returns: True if complete.

)doc";

static const char *__doc_separate_surface_overlapp =
R"doc(Separates two overlapping surfaces by contraction of surfaces boundary. 

Separates two overlapping surfaces iteratively by moving the vertices in the negative normal direction that is determined by the multiplication of a negative value and the shortest edge length corresponding to each vertex. This continues until the closest vertex is adjacent. The centeroids of the surfaces should be sufficiently
apart.

:param surf1: :class:`Surface` object.
:param surf2: :class:`Surface` object.
:param edge_movement: Multiplier of the smallest connected edge of a vertex, which equals the displacement of the vertices for each iteration. 
:param smoothing: Laplacian smoothing factor for each increment, see :func:`laplacian_smoothing`.
:param max_iter: The maximum number of iterations. 

:Returns: True if complete.

)doc";

static const char *__doc_separate_surface_overlapp_2 =
R"doc(Separates two overlapping surfaces outside a third surface.

Separates the points of two surfaces that are outside a third surface iteratively by moving the vertices in the negative normal direction that is determined by the multiplication of a negative value and the shortest edge length of each vertex. This continues until the closest vertex is adjacent/connected. The centeroids of the surfaces should be sufficiently apart. 

:param surf1: :class:`Surface` object.
:param surf2: :class:`Surface` object.
:param other: :class:`Surface` object, the algorithm is not applied to vertices inside this surface.
:param edge_movement: Multiplier of the smallest connected edge of a vertex, that gives the displacement of the vertices for each iteration. 
:param smoothing: Laplacian smoothing factor for each increment, see :func:`laplacian_smoothing`.
:param max_iter: The maximum number of iterations.
 
:Returns: True if complete.
    
)doc";

static const char *__doc_union_partially_overlapping_surfaces =
R"doc(Takes surface union of surfaces that partially overlapp each other.

Requires that the surfaces overlap. First, surface vertices inside the other surface are found. Then adjacent vertices within a threshold angle of the vertex normal are added to the vertex vector. The points corresponding to the vertex vector are iteratively moved in the vertex normal direction until there is sufficent overlapp between surfaces. 

:param surf1: :class:`Surface` object. 
:param surf2: :class:`Surface` object.
:param angle_in_degree: Threshold angle for vertex normal clustering, i.e. vertices with normal more than the threshold angle is not included in cluster.
:param adjustment: Multiplier of the smallest connected edge of a vertex, which equals the displacement of the vertices for each iteration. 
:param smoothing: Laplacian smoothing factor for each increment, see :func:`laplacian_smoothing`.
:param max_iter: The maximum number of iterations. 

:Returns: :class:`Surface` object.

)doc";

#if defined(__GNUG__)
#pragma GCC diagnostic pop
#endif

