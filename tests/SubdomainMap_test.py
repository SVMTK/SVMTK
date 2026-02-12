
import unittest
import SVMTK


class SubdomainMap_Test(unittest.TestCase):

    def test_add_subdomain(self):
        smap = SVMTK.SubdomainMap(2) 
        smap.add("01",2)
        
        a = smap.get_tags()
        self.assertEqual(a[0],0)         
        self.assertEqual(a[1],2) 
        smap.add("11",3)
        a = smap.get_tags()    
        self.assertEqual(a[2],3) 
        del smap

    def test_asterix_only(self):    
        smap = SVMTK.SubdomainMap(2) 
        smap.add("*",2)   
        bitmap = smap.get_map()
        self.assertEqual( bitmap['00'],0)             
        self.assertEqual( bitmap['01'],2)     
        self.assertEqual( bitmap['10'],2)     
        self.assertEqual( bitmap['11'],2)    
        del smap

    def test_dash_prefix(self):
        smap = SVMTK.SubdomainMap(3) 
        smap.add("-10",2)
        smap.add("-1",1)
        bitmap = smap.get_map()
        self.assertEqual( bitmap['010'],2)    
        self.assertEqual( bitmap['001'],1) 
        del smap
        
    def test_dash_suffix(self):  
        smap = SVMTK.SubdomainMap(3) 
        smap.add("1-",2)
        smap.add("01-",1)
        bitmap = smap.get_map()
        self.assertEqual( bitmap['100'],2)    
        self.assertEqual( bitmap['010'],1) 
        del smap
  
    def test_dash_and_astrix(self):     
        smap = SVMTK.SubdomainMap(3) 
        smap.add("*1-",2)
        bitmap = smap.get_map()
        self.assertEqual( bitmap['110'],2)    
        self.assertEqual( bitmap['010'],2)    
        del smap         
 
    def test_asterix_suffix(self):
        smap = SVMTK.SubdomainMap(3) 
        smap.add("100",2)
        smap.add("01*",3)
        bitmap = smap.get_map()
        self.assertEqual( bitmap['100'],2)    
        self.assertEqual( bitmap['010'],3)    
        self.assertEqual( bitmap['011'],3)    
        self.assertEqual( bitmap['000'],0)     
        del smap                    
     
    def test_asterix_init_error(self):            
        flag = False
        smap = SVMTK.SubdomainMap(0)
        try: 
            smap.add("*1",1)
        except SVMTK.InvalidArgumentError:
            flag = True
        self.assertTrue(flag)      
        del smap
   
    
    def test_erase_tag(self):
        smap = SVMTK.SubdomainMap(3) 
        smap.add("100",2)
        bitmap = smap.get_map()
        self.assertTrue('100' in bitmap )     
        smap.erase('100')
        bitmap = smap.get_map()
        self.assertFalse('100' in bitmap )    
        del smap


    def test_short_bitstring_error(self):    
        smap = SVMTK.SubdomainMap(3) 
        flag = False
        try: 
            smap.add("01",2)
        except SVMTK.InvalidArgumentError:
            flag = True
        self.assertTrue(flag)        
        del smap
     
            
    def test_long_bitstring_error(self):    
        smap = SVMTK.SubdomainMap(2) 
        flag = False
        try: 
            smap.add("00001",2)
        except SVMTK.InvalidArgumentError:
            flag = True
        self.assertTrue(flag)   
        del smap   
          
          
    def test_make_interfaces(self):
        s1 = SVMTK.Surface()
        s1.make_sphere(0.0,1.5,0.0,3.0,1.0)
        s2 = SVMTK.Surface()
        s2.make_sphere(0.0,-1.5,0.0,3.0,1.0)
        smap = SVMTK.SubdomainMap(2) 
        domain = SVMTK.Domain([s1,s2])
        domain.create_mesh(1.0)
        patches = domain.get_patches()
        interfaces = smap.make_interfaces(patches)
        self.assertEqual( len(interfaces), len(patches)) 
        smap.add_interface((1,0),9)  
        interfaces = smap.make_interfaces(patches)
        self.assertEqual( interfaces[(1,0)],9)
        del s1 
        del s2
        del smap
        
                  
             
    def test_add_interface(self):        
        smap = SVMTK.SubdomainMap(2) 
        smap.add_interface((1,0),2) 
        interfaces = smap.get_interfaces()
        print(interfaces)
        self.assertEqual(interfaces[(1,0)],2)  
        del smap    
                   

if __name__ == '__main__':
    unittest.main()




