How to fix mesh connectivity errors. 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


The connectivity errors can be divided into two groups, bad edge and bad vertex.

Bad edge is a boundary edge that is shared by more than 2 facets, which mostly occurs in narrow gaps that can not be resolved with the mesh resolution. This means that bad edges often occurs in 
bridges between narrow gaps in the volumetric mesh. 

(ADD IMAGE)

Bad vertex is a boundary vertex that is shared by two non-adjacent cells. This occurs with surface self-intersections and with a mismatch between surface and mesh resolution. (test Later)

(ADD IMAGE) 


We will start with how to fix bad edges. 
Removal of bad edges after the construction of the mesh can be volatile, so the strategi is to solve it by preprocessing the surfaces before constructing the mesh. This can be done by utilizing the function :func:`SVMTK.Surface.separate_narrow_gaps`, and you can follow other steps metioned in ( guide - How to repair surfaces). 

To fix bad vertices, we can use the function :func:`SVMTK.Surface.repair_self_intersections` or follow the less invasive method descriped in ( guide - How to repair surfaces)

Altough we have repaired the surface, we can still get bad edges or bad vertices if the mesh resolution is lower than the surface resolution. Therefore, we need to either fit the surface resolution to the mesh resolution or vice versa.


There is also problems related to thin but be aware that this method thus not take into account thin 

(Add image thin)



( Move to another guide)

We have a less invasive, but more expensive, option to remove self-intersection, which is shown in the following script.

(INLCUDE SCRIPT)



We have mentioned that selecting a mesh resolution should be comparable to the surface resolution, and we can use the function :func:`SVMTK.Surface.get_mesh_resolution` (not implemented) 



The surface resolution can be obtain with the function :func:`SVMTK.Surface.get_resolution` and can be altered by using :func:`SVMTK.Surface.isotropic_remeshing`. 
However, to simplify  

.. raw:: latex

    \newpage
