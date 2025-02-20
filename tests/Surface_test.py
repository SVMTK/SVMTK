import os
import unittest

import SVMTK
import pytest

tests_dir = os.path.dirname(__file__)


def ellipsoid_function( x, y, z):
  return x*x + 4.*y*y +4.*z*z-1.;

def torus_function(x,y,z):
     return ((x**2+y**2)**0.5 - 1)**2 + z**2 -4  
 
class Surface_Test(unittest.TestCase):

    def test_surface_io(self):
        surface =SVMTK.Surface()  
        surface.make_cube(0.,0.,0.,1.,1.,1.,1) 
        surface.save(f"{tests_dir}/Data/cube.off")
        surface.save(f"{tests_dir}/Data/cube.stl")
        surface= SVMTK.Surface(f"{tests_dir}/Data/cube.off")
        self.assertTrue(surface.num_vertices()==14 and surface.num_faces()==24 and surface.num_edges()==36) 
        surface= SVMTK.Surface(f"{tests_dir}/Data/cube.stl") 
        self.assertTrue(surface.num_vertices()==14 and surface.num_faces()==24 and surface.num_edges()==36)
        del surface

    def test_surfaces_shapes(self):
        surface =SVMTK.Surface()  
        surface.make_cube(0.,0.,0.,1.,1.,1.,1.) 
        self.assertTrue(surface.num_vertices()==14 and surface.num_faces()==24 and surface.num_edges()==36) 
        surface.clear() 
        self.assertEqual( surface.num_vertices(), 0) 
        surface.make_cone(0.,0.,0.,0,0.,2., 4.0,2.0, 8.37)
        self.assertTrue(surface.num_vertices()==8 and surface.num_faces()==12 and surface.num_edges()==18) 
        surface.clear() 
        self.assertEqual( surface.num_vertices(), 0) 
        surface.make_cone(0.,0.,0.,0,0.,2.,1.0,0.0,0.7) 
        self.assertTrue(surface.num_vertices()==25 and surface.num_faces()==46 and surface.num_edges()==69) 
        surface.clear() 
        self.assertEqual( surface.num_vertices(), 0) 
        surface.make_cylinder(0.,0.,0.,1.,1.0,.1,2.0,3.14) 
        self.assertTrue(surface.num_vertices()==10 and surface.num_faces()==16 and surface.num_edges()==24) 
        surface.clear() 
        self.assertEqual( surface.num_vertices(), 0) 
        del surface
    
    def test_surface_remeshing(self):
        surface =SVMTK.Surface()  
        surface.make_cube(0.,0.,0.,1.,1.,1.,1.) 
        self.assertTrue(surface.num_vertices()==14 and surface.num_faces()==24 and surface.num_edges()==36) 
        surface.isotropic_remeshing(1.0,1,1)
        self.assertTrue(surface.num_vertices()==14 and surface.num_faces()==24 and surface.num_edges()==36)
        del surface
        
    def boolean_operations(self):
        surface1 =SVMTK.Surface()  
        surface1.make_cube(0.,0.,0.,1.,1.,1.,1) 
        surface2 =SVMTK.Surface()  
        surface2.make_cube(0.5,0.5,0.5,1.5,1.5,1.5,1) 
        surface1.union(surface2) 
        self.assertTrue(surface.num_vertices()==20 and surface.num_faces()==36 and surface.num_edges()==54) 
        surface1.clear() 
        surface1.make_cube(0.,0.,0.,1.,1.,1.,1) 
        surface1.intersection(surface2) 
        self.assertTrue(surface.num_vertices()==8 and surface.num_faces()==12 and surface.num_edges()==18) 
        surface2.difference(surface1)        
        self.assertTrue(surface.num_vertices()==14 and surface.num_faces()==24 and surface.num_edges()==36) 
        del surface1
        del surface2
        
    def test_surface_meshing(self):
        surface = SVMTK.Surface();
        surface.make_sphere(0.0,0.0,0.0,1.0,1.0) 
        self.assertTrue(surface.num_vertices()>0 and surface.num_faces()>0 and surface.num_edges()>0)
        surface.clear()
        surface.implicit_surface(ellipsoid_function, 8 , 10)
        self.assertTrue(surface.num_vertices()>0 and surface.num_faces()>0 and surface.num_edges()>0) 
        del surface

    def test_span(self):
        surface =SVMTK.Surface()  
        surface.make_cube(0.,0.,0.,2.,1.,3.,1.)       
        self.assertEqual(surface.span(0),(0.,2))
        self.assertEqual(surface.span(1),(0.,1))
        self.assertEqual(surface.span(2),(0.,3))
        del surface

    def test_fill_holes(self):          
        mech_shark = SVMTK.Surface(f"{tests_dir}/Data/mech-holes-shark.off")
        closed, nb_holes = mech_shark.fill_holes()
        self.assertEqual(nb_holes,4)
        closed, nb_holes = mech_shark.fill_holes()
        self.assertEqual(nb_holes,0)
        del mech_shark 

    def test_triangulate_faces(self):
        surface = SVMTK.Surface(f"{tests_dir}/Data/P.off")
        self.assertTrue(surface.triangulate_faces())
        del surface

    def test_surface_clip(self):
        surface =SVMTK.Surface()
        surface.make_cube(-1.,-1.,-1.,1.,1.,1.,1)
        surface.clip(0,0,1.,0,True)
        self.assertAlmostEqual(surface.span(0)[0],-1.0,8)
        self.assertAlmostEqual(surface.span(0)[1],1.0,8)
        self.assertAlmostEqual(surface.span(1)[0],-1.0,8)
        self.assertAlmostEqual(surface.span(1)[1],1.0,8)
        self.assertAlmostEqual(surface.span(2)[0],-1.,8)
        self.assertAlmostEqual(surface.span(2)[1],0.0,8)
        del surface

    def test_adjust_boundary(self):
        surface =SVMTK.Surface()   
        surface.make_cube(-1.,-1.,-1.,1.,1.,1.,1) 
        surface.adjust_boundary(-0.1)
        self.assertTrue(surface.span(0)[1]-surface.span(0)[0] < 2. ) 
        surface.clear() 
        surface.make_cube(-1.,-1.,-1.,1.,1.,1.,1) 
        surface.adjust_boundary(0.1)
        self.assertTrue(surface.span(0)[1]-surface.span(0)[0] > 2. )
        del surface
    
    def test_smoothing(self):
        surface =SVMTK.Surface()   
        surface.make_cube(-1.,-1.,-1.,1.,1.,1.,1) 
        self.assertAlmostEqual(surface.volume(),8)
        surface.smooth_laplacian(0.8,1) 
        self.assertTrue(surface.volume()<7)
        surface.make_cube(-1.,-1.,-1.,1.,1.,1.,1) 
        surface.smooth_taubin(1) 
        self.assertTrue(surface.volume()>7)
        del surface
        

    def test_mean_curvature_flow(self):
        surface =SVMTK.Surface()
        surface.implicit_surface(torus_function, 6,30) 
        l1 =surface.mean_curvature_flow()
        self.assertTrue( l1[0]==l1[-1])
        del surface
        
    def test_shortest_surface_path(self):      
        surface =SVMTK.Surface()   
        surface.make_cube(-1.,-1.,-1.,1.,1.,1.,1)         
        a = surface.get_shortest_surface_path(SVMTK.Point_3(-1,-1,1), SVMTK.Point_3(1,-1,1))
        self.assertTrue( a[-1]==SVMTK.Point_3(-1,-1,1) and a[0]==SVMTK.Point_3(1,-1,1))
        del surface

    def test_collapse_and_split_edges(self): 
        surface =SVMTK.Surface()  
        surface.make_cube(0.,0.,0.,1.,1.,1.,1) 
        self.assertTrue(surface.num_vertices()==14 and surface.num_faces()==24 and surface.num_edges()==36) 
        surface.split_edges(0.5)
        self.assertEqual(surface.num_edges(),144) 
        surface.collapse_edges()
        self.assertEqual(surface.num_edges(),18) 
        surface.split_edges(0.5)
        surface.collapse_edges(1.0)
        self.assertEqual(surface.num_edges(),6) 
        del surface

    def test_extension(self):
        surface=SVMTK.Surface()   
        surface.make_cube(-1.,-1.,-1.,1.,1.,1.,0.5) 
        cylinder = surface.extension(SVMTK.Point_3(0,0,2),0.5,0.5,0.25,True)
        self.assertTrue(cylinder.num_edges()>0) 
        del surface

    def test_separate_narrow_gaps(self):
        surface=SVMTK.Surface(f"{tests_dir}/Data/narrow_gap.off")
        result = surface.separate_narrow_gaps(-0.5,0.5)
        self.assertTrue( result[0]) 
        self.assertEqual(result[1],0) 
        del surface

    def test_separate_close_vertices(self):
        surface =SVMTK.Surface()  
        surface.set_proximity_ratio(1.0)
        surface.make_cube(0.,0.,0.,0.1 , 1., 1.,.01) 
        result = surface.separate_close_vertices(0.5)
        self.assertTrue(result[0]) 
        self.assertEqual(result[1],0) 
        del surface
         
    def test_enclose(self):
        surface1 =SVMTK.Surface()  
        surface1.make_cube(0.,0.,0.,2.,2.,2.,2.0) 
        surface2 =SVMTK.Surface()  
        surface2.make_cube(-0.5,-0.5,-0.5,2.5,2.5,2.5,2.) 
        a =  surface1.enclose(surface2,0.8,0.1) 
        self.assertTrue( a[0])
        self.assertEqual( a[1],0) 
        del surface1
        del surface2
        
    def test_embed(self):
        surface1 =SVMTK.Surface() 
        surface1.make_cube(0.,0.,0.,2.,2.,2.,2) 
        surface2 =SVMTK.Surface()  
        surface2.make_cube(-0.5,-0.5,-0.5,2.5,2.5,2.5,2.) 
        a =  surface2.embed(surface1,-0.8)
        self.assertTrue( a[0])
        self.assertEqual( a[1],0) 
        del surface1
        del surface2
        
    def test_separate_surface(self):
        surface1 =SVMTK.Surface()  
        surface1.make_cube(0.,0.,0.,2.,2.,2.,0.5)         
        surface2 =SVMTK.Surface()  
        surface2.make_cube(0.05,0.05,0.05,1.95,1.95,1.95,0.5) 
        a = surface2.separate(surface1, 1.0, 0.5 ,600)    
        self.assertTrue( a[0])
        self.assertEqual( a[1],0)
        del surface1
        del surface2
        
    def test_reconstruction(self):
        surface =SVMTK.Surface()  
        surface.make_cube(0.,0.,0.,10.,10.,10.,1.) 
        surface.reconstruct(20,1.0,1.0)
        surface.collapse_edges() 
        self.assertTrue(surface.num_vertices()>0 and surface.num_faces()>0 and surface.num_edges()>0)
        del surface

    def test_convex_hull(self):
        surface1=SVMTK.Surface()   
        surface1.make_cube(-1.,-1.,-1.,1.,1.,1.,1) 
        surface2 =surface1.convex_hull()
        self.assertEqual(surface2.num_vertices(),8)
        del surface1
        del surface2




if __name__ == '__main__':
    import os
    unittest.main()
    os.remove('tests/Data/cube.stl')
    os.remove('tests/Data/cube.off')
