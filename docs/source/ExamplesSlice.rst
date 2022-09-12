:orphan:

.. _ExampleSlicedoc:

Slice Examples
=================================

The SVMTK Slice class can be used to create a 2D mesh, and we will here show some examples on how to utilize the different functionalities.



How to make a 2D mesh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The first example is on how to construct a 2D mesh given some constraints, i.e. boundaries. 

.. literalinclude:: ../../examples/Slice/example_make_2d_mesh.py
   :lines: 1-22
   
.. figure:: ../Images/Hexagon.png
   :width: 40%
   :alt: Hexagon slice
   :align: center

   .. 
   
   Hexagon slice mesh.

In some cases, contradicting constraints can cause the mesh construction to crash, as the mesh criteria fails to create a mesh with both constraints. This would happen in the next example, however if we slightly change the hexagon radius, it will no longer crash.

.. literalinclude:: ../../examples/Slice/example_make_2d_mesh.py
   :lines: 24-


.. figure:: ../Images/Flower.png
   :width: 40%
   :alt: Hexagon slice with petals.
   :align: center
   
   .. 
   
   Hexagon slice mesh with petals.


Combining Slices
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example will show the method to combine slice constraints from different slices.

.. literalinclude:: ../../examples/Slice/example_combine_slices.py

   
.. figure:: ../Images/Combine_slices.png
   :width: 40%
   :alt: Slice with 3 subdomains
   :align: center

   .. 
   
   Combining a circle and a square slice.

Add Subdomains to a slice mesh.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this example, we will show how to add and also remove subdomains from a slice mesh. The methodology is similar to :class:`SVMTK.Surface`
and it works with :class:`SVMTK.SubdomainMap` 

.. literalinclude:: ../../examples/Slice/example_add_subdomain_tags.py



.. figure:: ../Images/Slice_3subdomains.png 
   :width: 40%
   :alt: Slice with 3 subdomains
   :align: center
   
   ..
   
   Adding 3 subdomains to the mesh.
   
   
.. figure:: ../Images/Slice_2subdomains.png 
   :width: 40%
   :alt: Union boolean operation
   :align: center
   
   ..
   
   Removing one subdomains from the mesh.
   

Keeping largest component 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Slicing complex surfaces can result in isolated clusters of faces. In this example, we will show the function :func:`SVMTK.Slice.keep_largest_connected_component`. It is worth mentioning that this method require the function :func:`SVMTK.Slice.add_surface_domains` to remove helping faces.

.. literalinclude:: ../../examples/Slice/example_keep_largest_connected_component.py


.. figure:: ../Images/Unfethered.png
   :width: 80%
   :alt: Helping faces
   :align: center
   
   .. 
   
   Slice mesh with helping faces.


.. figure:: ../Images/keep_all_components.png 
   :width: 80%
   :alt: keep all components
   :align: center
   
   .. 
   
   Slice mesh without helping faces.

.. figure:: ../Images/keep_largest_component.png 
   :width: 80%
   :alt: keep largest component
   :align: center
   
   .. 
   
   Slice mesh of the largest connected component.


   
Export as surface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The z-coordinate in slice meshes will be zero. (Likely to be removed, need additional utility) 


.. figure:: ../Images/export_as_surface_2.png   
   :width: 50%
   :alt: Smoothing Laplacian surface
   :align: center 

   .. 
   
   The first and second slice mesh respectively marked as yellow and red.


.. raw:: latex

    \newpage

