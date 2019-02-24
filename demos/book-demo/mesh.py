#!/usr/bin/python

if __name__ =='__main__':
        import SVMTK as svm
        import argparse
        import os

        parser = argparse.ArgumentParser(description='Creates a volume mesh of an input surface' )

        parser.add_argument("--i"   ,type=str, help = "The path to the input surface.",required=True)
        parser.add_argument("--res"   ,type=float, default=32.0    ,help = "The resoltuion of the mesh. Defined as cell_size = bounding_radius/resolution")
        parser.add_argument("--o"   ,type=str, help = "The path to the output mesh.")
        #parser.add_argument("--proc"  ,action='store_true',help = "If included, the script will use the {}-{}-proc.stl as input.".format(hemi,surface))
        Z = parser.parse_args()
 
        surf  = svm.Surface(Z.i)

        domain = svm.Domain(surf)
        domain.create_mesh(Z.res)
 
        domain.save(Z.o)



