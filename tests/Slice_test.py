import os
import unittest

import SVMTK

tests_dir = os.path.dirname(__file__)


class Slice_Test(unittest.TestCase):


    def test_unpack_mesh(self): 
        slc = SVMTK.Slice();
        slc.add_constraint([SVMTK.Point_2(0,0),SVMTK.Point_2(0,1),SVMTK.Point_2(1,1), SVMTK.Point_2(1,0),SVMTK.Point_2(0,0) ] )
        slc.create_mesh(8.0)
        points    = slc.get_points()
        triangles = slc.get_cells() 
        edges     = slc.get_facets()  
        self.assertEqual(points.shape,    (365, 2))        
        self.assertEqual(triangles.shape, (664, 3))
        self.assertEqual(edges.shape,     ( 64, 2))

    def test_slice_add_constrains(self):
        slice1 = SVMTK.Slice()
        slice1.add_constraint([SVMTK.Point_2(0,0),SVMTK.Point_2(1,1)]) 
        self.assertTrue( slice1.num_constraints(),1)
        slice2 = SVMTK.Slice(SVMTK.Plane_3(0.,0.,1.,0.) )
        slice2.add_constraint([SVMTK.Point_2(0.5,0.5),SVMTK.Point_2(0.5,-0.5)]) 
        slice3 = slice2
        slice3.add_constraints(slice1)
        self.assertTrue( slice3.num_constraints(),2)

    def test_slice_meshing(self):
        surface = SVMTK.Surface() 
        surface.make_cube(0.,0.,0.,1.,1.,1.,1) 
        slice_ = surface.get_slice(1,1,1,-0.5 )  
        self.assertTrue( slice_.num_constraints() > 0)
        slice_.create_mesh(5.) 
        self.assertTrue( slice_.num_cells() > 0) 
    

    def test_slice_subdomains(self):
        slice_ = SVMTK.Slice(SVMTK.Plane_3(0,0,1,0))
        surface1 = SVMTK.Surface() 
        surface1.make_cube(-1.,-1.,-1.,1.,1.,1.,2.) 
        surface2 = SVMTK.Surface() 
        surface2.make_cube(-2.,-2.,-2.,2.,2.,2.,4.)
        
        sf= SVMTK.SubdomainMap(0) 
        sf.add("11",2)
        sf.add("11",2)
        slice_ = SVMTK.Slice(SVMTK.Plane_3(0,0,1,0))
        slice_.slice_surfaces([ surface1 , surface2] )  

        slice_.create_mesh(5.)
        slice_.add_surface_domains([ surface1 , surface2] )
        self.assertEqual(slice_.num_subdomains(),2)
        slice_.add_surface_domains([ surface1 , surface2],sf)
        self.assertEqual(slice_.num_subdomains(),1)
        slice_.remove_subdomain(2)
        try:
          self.assertEqual(slice_.num_subdomains(),0)
        except SVMTK.EmptyMeshError:
          self.assertTrue(True) 

    def test_connected_components(self): 
        surface1 = SVMTK.Surface() 
        surface1.make_cube(0.,0.,0.,1.,1.,1.,0.1)
        surface2 = SVMTK.Surface() 
        surface2.make_cube(1.5,0,0,2.,1.,1.,0.1)
        slice_ = SVMTK.Slice(SVMTK.Plane_3(0,0,1,-0.5))
        slice_.slice_surfaces([surface1,surface2])  
        slice_.create_mesh(6.)
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


    def test_constraint_tags(self): 
        slc = SVMTK.Slice() 
        slc = SVMTK.Slice()
        slc.add_constraint([SVMTK.Point_2(0,0),SVMTK.Point_2(0,1),SVMTK.Point_2(1,1), SVMTK.Point_2(1,0),SVMTK.Point_2(0,0) ] ) 
        slc.create_mesh(1.0) 
        self.assertEqual(slc.get_facet_tags().shape[0],32)
    

    def test_slice_mesh(self):
        domain = SVMTK.Domain(f"{tests_dir}/Data/cube.mesh") 
        slc = SVMTK.Slice(0,0,1,-0.5) 
        slc.slice_mesh(domain) 
        slc.create_mesh(5.0) 
        self.assertEqual( slc.get_cells().shape ,(430, 3) )

if __name__ == '__main__':
    unittest.main()
    import os

    os.remove(f"{tests_dir}/Data/slice.vtu")
    os.remove(f"{tests_dir}/Data/slice.stl")
    os.remove(f"{tests_dir}/Data/slice.off")
    os.remove(f"{tests_dir}/Data/slice.mesh")
