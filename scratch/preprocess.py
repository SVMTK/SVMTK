import sys

from dolfin_utils.meshconvert import xml_writer 
import re 


def run(name ): 
  """ 
  convert between amira and fenics implemented as state machine 
  0 = look for starters 
  1 = vertex 
  2 = cells 
  3 = labels
  """
  state = 0 
  ifile = open(name+".mesh")
  ofile = open(name+".xml", "w")
  ofile_sub = open(name+"_sub.xml", "w")
  
  vertices =[] 
  cells = []
  labels = []
  coordinates =[]
  index_vertices = 0
  index_cells = 0
  while 1:
    line = ifile.readline()
    if line.startswith("Vertices"): 
      state = 1
    elif line.startswith("Triangles"): 
      state = 2  
    elif line.startswith("Tetrahedra"):
      state = 3
    elif line.startswith("End"):
      break
    
    elif state == 1:   
    # example: 
    #-37.168381 50.233715 -12.345407
      #print line.strip().split()
      if len (line.strip().split()) > 1:
	      index_vertices += 1
	      x, y, z, label = line.strip().split()
              x = float(x)
	      y = float(y)
              z = float(z)
	      vertices.append((index_vertices,x,y,z))
              coordinates.append((x,y,z))

    elif state == 3:   
    # example: 
#    *tetra4(1,1,38366,38367,9253,43441)
      if len (line.strip().split()) > 1:
	      index_cells += 1
	      v1, v2, v3, v4, label= line.strip().split() 
	      cells.append((index_cells,v1,v2,v3,v4))
	      labels.append(label)

  from collections import Counter

  print len(coordinates) == len(set(coordinates))
  print "out of loop"
  xml_writer.write_header_mesh(ofile, "tetrahedron", 3) 

  xml_writer.write_header_vertices(ofile, len(vertices))
  counter = 0 
  for v in vertices: 
    if not counter == int(v[0])-1: print "Warning! Numbering assumption invalidated", counter, v[0]  
    xml_writer.write_vertex(ofile, int(v[0])-1, float(v[1]), float(v[2]), float(v[3])) 
    counter += 1 
  xml_writer.write_footer_vertices(ofile)

  xml_writer.write_header_cells(ofile, len(cells))
  counter=0
  for c in cells: 
    if not counter == int(c[0])-1: print "Warning! Numbering assumption invalidated", counter, c[0]  
    xml_writer.write_cell_tetrahedron(ofile, int(c[0])-1, int(c[1])-1, int(c[2])-1, int(c[3])-1, int(c[4])-1) 
    counter += 1 
  xml_writer.write_footer_cells(ofile)


  header = ("""<?xml version="1.0"?>
<dolfin xmlns:dolfin="http://fenicsproject.org">
  <mesh_function>
    <mesh_value_collection name="f" type="uint" dim=\"%d" size=\"%d">\n""") % (3, len(cells))
  
  ofile_sub.write(header)
  
  counter=0
  for c in cells: 
    if not counter == int(c[0])-1: print "Warning! Numbering assumption invalidated", counter, c[0]
    xml_writer.write_entity_meshvaluecollection(ofile_sub, 1, int(c[0])-1, int(labels[counter]), local_entity=0) 
    counter += 1
    
  xml_writer.write_footer_meshvaluecollection(ofile_sub)
  xml_writer.write_footer_meshfunction(ofile_sub)
  xml_writer.write_footer_mesh(ofile)
 
#   not handled: 
#    *component(1,"?",1,1)
#    *tria3(1,1,9253,38366,38367)
#
#    *component(2,"Inside",1,2)
#    *tetra4(1,1,38366,38367,9253,43441)
#      
#    if line.startswith("

if __name__ == "__main__": 
  import sys
  print sys.argv[1]
  run(sys.argv[1])

