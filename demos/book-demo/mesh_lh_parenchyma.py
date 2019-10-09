#!/usr/bin/python


import SVMTK as svm 



if __name__ =='__main__':

        import argparse
        import os

        parser = argparse.ArgumentParser(description='Creates a volume mesh of a hemisphere parenchyma, given pial, white and ventricluar surfaces. It is assumed that the pial encloses both ventriclea and white surfaces.' )

        parser.add_argument("--pial"       ,type=str, help = "The path to the pial surface", required=True)
        parser.add_argument("--white"      ,type=str, help = "The path to the white surface",required=True)
        parser.add_argument("--ventricle"  ,type=str, help = "The path to the ventricle surface",required=True)
        parser.add_argument("--res"        ,type=float, default=32.0    ,help = "The resoltuion of the mesh. Defined as the upperbound of the cell_size, equal to bounding_radius/resolution. The default value is 32.")
        parser.add_argument("--keep"       ,action='store_true',help = "Include to mark and keep ventricle surface.")

        parser.add_argument("--surf"       ,type=str, help = "The path to the output mesh boundary.")
        parser.add_argument("--o"          ,type=str, help = "The path to the output mesh.",required=True)


        Z = parser.parse_args()
 
        pial  = svm.Surface(Z.pial)             # Loads the pial surface into a Surface Obejct.
        white = svm.Surface(Z.white)            # Loads the white surface into a Surface Obejct.
        vent  = svm.Surface(Z.ventricle)        # Loads the ventricle surface into a Surface Obejct.


        surfaces = [pial,white ,vent]
        smap      = svm.SubdomainMap()          # Creates an empty SubDomainMap Obejct.
        smap.add("100",1)                       # Marks all cells that are only inside the pial surface.
        smap.add("110",2)                       # Marks all cells inside both the pial and white surface.
       
        
        smap.add("111",3)                    # Marks all cells inisde the pial, white and ventricle surface.

 
        domain   = svm.Domain(surfaces,smap)    # Creates a domain with specified subdomains
       
        domain.create_mesh(Z.res)               # Generates the volume mesh of the loaded surface.
  
        domain.remove_subdomain(3)              # removes cells with subdomain 3 tag 
        domain.save(Z.o)                        # Save the mesh at the 

        if Z.surf is not None:
           surf = domain.get_boundary()         # Extracts all boundary facets to a SVM Surface Object
           surf.save(Z.surf)                         # Saves the Surface Object containing the boundary
