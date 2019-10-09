#!/usr/bin/python

if __name__ =='__main__':
        import SVMTK as svm
        import argparse
        import os

        parser = argparse.ArgumentParser(description='Creates a volume mesh of an input surface' )

        parser.add_argument("--i"   ,type=str, help = "The path to the input surface. Both .off and .stl can be loaded",required=True)
        parser.add_argument("--res"   ,type=float, default=32.0    ,help ="The resoltuion of the mesh. Defined as the upperbound of the cell_size, equal to bounding_radius/resolution. The default value is 32.")
        parser.add_argument("--o"   ,type=str, help = "The path to the output mesh. Requires extension .mesh")
        #parser.add_argument("--proc"  ,action='store_true',help = "If included, the script will use the {}-{}-proc.stl as input.".format(hemi,surface))
        Z = parser.parse_args()
 
        surf  = svm.Surface(Z.i)            # Loads the surface into the Surface obejct

        domain = svm.Domain(surf)           # Loads the Surface Object into the Domain object
        
        domain.create_mesh(Z.res)           # Generates the volume mesh of the loaded surface.
 
        domain.save(Z.o)                    # Save the mesh at the 


 
