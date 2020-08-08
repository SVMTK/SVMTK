import unittest
import SVMTK


class Slice_Test(unittest.TestCase):


    def test_slice_add_constrains(self):
        slice1 = SVMTK.Slice()
        slice1.add_constraint([SVMTK.Point_2(0,0),SVMTK.Point_2(1,1)]) 
        self.assertTrue( slice1.number_of_constraints(),1)
        slice2 = SVMTK.Slice(SVMTK.Plane_3(0.,0.,1.,0.) )
        slice2.add_constraint([SVMTK.Point_2(0.5,0.5),SVMTK.Point_2(0.5,-0.5)]) 
        slice3 = SVMTK.Slice(slice2)
        slice3.add_constraints(slice1)
        self.assertTrue( slice3.number_of_constraints(),2)

    def test_slice_meshing(self):
        surface = SVMTK.Surface() 
        surface.make_cube(0.,0.,0.,1.,1.,1.,1) 
        slice_ = surface.slice(1,1,1,-0.5 )  
        self.assertTrue( slice_.number_of_constraints() > 0)
        slice_.create_mesh(1.) 
        self.assertTrue( slice_.number_of_faces() > 0) 
        slice_.save("tests/Data/slice.vtu")
        slice_.save("tests/Data/slice.stl")
        slice_.save("tests/Data/slice.off")
        slice_.save("tests/Data/slice.mesh")


    def test_slice_subdomains(self):
        slice_ = SVMTK.Slice(SVMTK.Plane_3(0,0,1,0))
        surface1 = SVMTK.Surface() 
        surface1.make_cube(-1.,-1.,-1.,1.,1.,1.,1) 
        surface2 = SVMTK.Surface() 
        surface2.make_cube(-2.,-2.,-2.,2.,2.,2.,1)


        sf= SVMTK.SubdomainMap() 
        sf.add("11",2)
        sf.add("11",2)
        slice_ = SVMTK.Slice(SVMTK.Plane_3(0,0,1,0))
        slice_.slice_surfaces([ surface1 , surface2] )  

        slice_.create_mesh(1.)
        slice_.add_surface_domains([ surface1 , surface2] )
        self.assertEqual(slice_.number_of_subdomains(),2)
        slice_.add_surface_domains([ surface1 , surface2],sf)
        self.assertEqual(slice_.number_of_subdomains(),1)
        slice_.remove_subdomains(2)
        self.assertEqual(slice_.number_of_subdomains(),0)

    def test_connected_components(self):
        surface1 = SVMTK.Surface() 
        surface1.make_cube(0.,0.,0.,1.,1.,1.,10)
        surface2 = SVMTK.Surface() 
        surface2.make_cube(1.5,0,0,2.,1.,1.,10)
        slice_ = SVMTK.Slice(SVMTK.Plane_3(0,0,1,-0.5))
        slice_.slice_surfaces([ surface1,surface2])  
        slice_.create_mesh(1.)
        slice_.add_surface_domains([ surface1,surface2])
        self.assertEqual(slice_.connected_components(), 2) 
        slice_.keep_largest_connected_component() 
        self.assertEqual(slice_.connected_components(), 1) 


    def test_simplify(self):
        slice1 = SVMTK.Slice()
        slice1.add_constraint([SVMTK.Point_2(0,0),SVMTK.Point_2(1,1),SVMTK.Point_2(2,1)] ) 
        slice1.simplify(1.0) 
        constraints = slice1.get_constraints()
        self.assertEqual(len(constraints[0]), 3) 
 
    #def tearDown(self):
    #    import os
    #    os.remove("./Data/slice.vtu")
    #    os.remove("./Data/slice.stl")
    #    os.remove("./Data/slice.off")
    #    os.remove("./Data/slice.mesh")


if __name__ == '__main__':
    unittest.main()

