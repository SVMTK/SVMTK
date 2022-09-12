""" Create 3D mesh with two conforming edge helixes inside """ 

import SVMTK as svm
import numpy as np
from pathlib import Path

if __name__ == "__main__":
   print("Start ",__file__)       
   surf = svm.Surface();
   outdir = Path("results")
   outdir.mkdir(exist_ok=True)
   # Lower center coordinates 
   x0,y0,z0 = 0,0,0
   # Upper center coordinates
   x1,y1,z1 = 0,0,10
   # Cylinder radius 
   r = 5 
   # Edge length 
   edge_length = 1.
   # Create cylinder
   surf.make_cylinder(x0,y0,z0,x1,y1,z1,r,edge_length)
   # Helix radius
   R=3.0
   helix1=[]
   helix2=[]
   for i in np.linspace(0.5,9.5,40) : 
       helix1.append(svm.Point_3( R*np.cos(np.pi*i*0.5), R*np.sin( np.pi*i*0.5) ,i) )
       helix2.append(svm.Point_3(-R*np.cos(np.pi*i*0.5), -R*np.sin( np.pi*i*0.5),i) ) 
       
   maker = svm.Domain(surf)
   maker.add_feature(helix1)
   maker.add_feature(helix2)
   # Add the sharp edges of the cylinder surface.
   maker.add_sharp_border_edges(surf,70)
   maker.create_mesh(22.)
   maker.boundary_segmentations()
   maker.save(str(outdir/"Helixfeature.mesh"))
   print("Finish ",__file__)       

