import unittest
import SVMTK


class Slice_Test(unittest.TestCase):
    def test_slice_init(self):
        slice_ = SVMTK.Slice()
        slice_ = SVMTK.Slice(SVMTK.Plane_3(0.,0.,1.,0.) )
        slice_ = SVMTK.Slice(slice_)


    def test_slice_meshing(self):
        surface = SVMTK.Surface() 
        surface.make_cube(0.,0.,0.,1.,1.,1.,1) 
        slice_ = surface.slice(1,1,1,-0.5 )  

        self.assertTrue( slice_.number_of_constraints() > 0)
        slice_.create_mesh(1.) 
        self.assertTrue( slice_.number_of_faces() > 0) 


    def test_slice_surfaces(self):
        slice_ = SVMTK.Slice(SVMTK.Plane_3(0,0,1,0))
        slice_.slice_surfaces([SVMTK.Surface(),SVMTK.Surface()] )  

        slice_ = SVMTK.Slice(SVMTK.Plane_3(0,0,1,0))
        surface_1 = SVMTK.Surface() 
        surface_1.make_cube(-1.,-1.,-1.,1.,1.,1.,1) 
        surface_2 = SVMTK.Surface() 
        surface_2.make_cube(-2.,-2.,-2.,2.,2.,2.,1)
        sf= SVMTK.SubdomainMap() 

        slice_ = SVMTK.Slice(SVMTK.Plane_3(0,0,1,0))
        slice_.slice_surfaces([ surface_1 , surface_2] )  





if __name__ == '__main__':
    unittest.main()

