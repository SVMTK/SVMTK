""" Creating 2D mesh using constraints """  

import SVMTK as svm
import numpy as np
from pathlib import Path


if __name__ == "__main__":
   print("Start ",__file__)   
   outdir = Path("results")
   outdir.mkdir(exist_ok=True) 
   slc = svm.Slice()
   # Hexagon radius 
   r = 8   
   # Edges of a hexagon
   hexagon = [svm.Point_2( r*np.cos(rad),r*np.sin(rad)) for rad in np.deg2rad([0,60,120,180,240,300])]  
   # Closing the constraints  
   hexagon.append(hexagon[0])
   # Adding constraints to the slice 
   slc.add_constraint(hexagon)
   # Creating a mesh
   slc.create_mesh(32)
   slc.save("hexagon.vtu")

   # Slightly altering the hexagon radius.
   r = 7.99
   for n in range(0,6):
       rad = n*np.pi/3 # TODO SIMPLIFY
       petal = [svm.Point_2( r*( np.cos(rad) + np.cos(theta)) , r*(np.sin(rad) + np.sin(theta)) ) for theta in np.linspace( rad + 2*np.pi/3, rad + 4*np.pi/3 ,100)  ]
       slc.add_constraint(petal)                    

   # Creating a new mesh will clear the old mesh. 
   slc.create_mesh(32)
   slc.save("flower.vtu")
   print("Finish ",__file__)  
