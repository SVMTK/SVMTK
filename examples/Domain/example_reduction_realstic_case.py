""" Create a mesh with multiple surfaces """ 

import SVMTK as svm
import numpy as np
from pathlib import Path

if __name__ == "__main__":
   print("Start ",__file__)       
   surf = svm.Surface();
   outdir = Path("results")
   outdir.mkdir(exist_ok=True)


   pial = svm.Surface("../Data/lh-pial.stl") 
   
   convex_hull = pial.convex_hull()
   

   convex_hull.isotropic_remeshing(1.0,5,False) 
   
   convex_hull.adjust_boundary(2.0) 
   
   
   convex_hull.isotropic_remeshing(1.0,5,False) 
   
   convex_hull.save("pial_convex_hull.stl")
   exit()
   smap = svm.SubdomainMap()

   smap = svm.SubdomainMap()
   # Adding subdomains 
   smap.add("10",1)
   smap.add("11",2)   
   

   maker = svm.Domain([convex_hull,pial],smap)
   maker.create_mesh(32)


   maker.subdomain_reduction("sas_depth.dat",1)

   print(maker.number_of_facets() )
   

   
   
   
   
   
   maker.save("reduced_sas.mesh")
   print("Finish ",__file__)       

