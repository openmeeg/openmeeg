Each directory contains the files needed for describing the head model geometry, and a forward computation, i.e the sensors locations and the sources (dipoles..).

==========
 Models 
==========

Nested models:
^^^^^^^^^^^^^^
Head1, Head2, Head3 are three nested-sphere models, consisting of three spheres with radii 0.87, 0.92 and 1.

Model Head1 has 42 vertices per surface
Model Head2 has 162 vertices per surface
Model Head3 has 642 vertices per surface


Non-Nested models:
^^^^^^^^^^^^^^^^^^
- HeadNNa1, HeadNNa2, HeadNNa3, 
- HeadNNb1, HeadNNb2, HeadNNb3,
- HeadNNc1, HeadNNc2, HeadNNc3 are Non-Nested models.

- HeadNNaX are (rotated) HeadX, where the inner sphere is split in two hemispheres: north & south. ( with X={1,2,3})
- HeadNNcX are versions of HeadNNaX with a different (finer) discretization close to the two hemispheres borders: north & south. ( with X={1,2,3})
         HeadNNc1 and HeadNNc3 were generated using the CGAL tool provided, while model HeadNNc2 is one of the head models used in the article of Kybic et al [1]_ ( with X={1,2,3})

.. [1] Jan Kybic, Maureen Clerc, Olivier Faugeras, Renaud Keriven, Théo Papadopoulo. Generalized head models for MEG/EEG: boundary element method beyond nested volumes, Physics in Medicine and Biology, vol. 51, pp. 1333--1346, March 2006.

- HeadNNbX is HeadX, where in the inner sphere were added two little spheres of radius 0.3 with center at (0., +/-0.4, 0.): spherenorth & spheresouth. ( with X={1,2,3})

- Head0 and HeadNNa0 are test cases for non-regression testing.

================================
 Writting geom and cond files 
================================

\*.geom files:
^^^^^^^^^^^^^^

OpenMEEG definitions:
---------------------
"Mesh": a mesh is a collection of vertices and triangles, all connected.
"Interface": an interface is a closed mesh.
"Half-space": an interface then defines a Inside and Outside. An Half-space is thus a space INSIDE or OUTSIDE an interface.
"Domain": a domain is an intersection of half-spaces.


Namings:
--------
For ease of explanation, we consider the naming as follow:
- meshes: are named in lower case
- interfaces: in CamelCase
- domains: in UPPERCASE

but this convention is up to you.

A *.geom file* contains:

1. A header section: 
--------------------

Version 1.1 and later date from after the introduction of non-nested geometries::

  # Domain Description 1.1                             


2. An (optional) Meshes section: 
--------------------------------

It indicates which files contains the geometry::

  Meshes 3                                             
                                                      
  Mesh north: "north.tri"                              
  Mesh south: "south.tri"                              
  Mesh skull: "skull.tri"                              
  Mesh scalp: "scalp.tri"                              


Supported mesh formats are tri, off, bnd, mesh, vtk (provided openmeeg was compiled with USE_VTK=True).

Note that the individual meshes can either be named or not. If no name is provided, they will be automatically named '1', '2', ... following the order.

2. An (optional) MeshFile section: 
----------------------------------
In case you did not load the meshes through the Meshes section, and if VTK is enabled, the best format to use is VTK/vtp (at least for non-nested geometries)::

  MeshFile Head.vtp

In this format, all meshes are in a single VTK/vtp file, where all polygons (triangles) have a string data attached
indicating the name of the mesh it belongs to. (These files can easily be opened with Paraview www.paraview.org, select some triangles-> Cell Label-> check Visible, see HeadNNa1/HeadNNa1.png for an example of visualization.)


3. An Interfaces section:
-------------------------
Starting with the keyword "Interfaces", followed by the number of interfaces.

In case no MeshFile section was found, and no Meshes section neither, you have to use the format in the example below::

  Interfaces 3
  
  Interface: "cortex.1.tri"
  Interface: "skull.1.tri"
  Interface: "scalp.1.tri"
  

In the example below, meshes are re-oriented with a plus or minus sign, so that the overall interface (composed by several meshes or only one) is consistently oriented (OpenMEEG will complain if this is not the case)::
 
  Interfaces 5                               
                                            
  Interface North: +north +cut               
  Interface South: +south -cut               
  Interface Cortex: north south              
  Interface Skull: skull                     
  Interface Scalp: scalp                     

Note that the interfaces can either be named or not. If no name is provided, they will be automatically named '1', '2', ... following its order in the list.



4. A Domains section:
---------------------
The definition of the domains uses the previously defined interfaces::

  Domains 5                                  
                                           
  Domain NORTH: -North                       
  Domain SOUTH: -South                       
  Domain SKULL: -Skull +Cortex               
  Domain SCALP: -Scalp +Skull                
  Domain AIR: +Scalp                         

A Domain is defined as being OUTSIDE or INSIDE certain interfaces.


\*.cond files:
^^^^^^^^^^^^^^

A *.cond file* contains:

1. A header section: 
--------------------
::

  # Properties Description 1.0 (Conductivities) 

2. Conductivity values: 
-----------------------
::
 
  AIR         0.0                            
  NORTH       1                              
  SOUTH       1                              
  SKULL       0.0125                         
  SCALP       1                              
                                                
Each domain name is followed by its conductivity value.

============================================
 Example for generating meshes and vtp files 
============================================

Using the tools:
^^^^^^^^^^^^^^^^
CGAL_  can generate surfacic meshes out of implicit functions, 3D image levelsets,...
We here show how models such as *HeadNNc1* and *HeadNNc3* were generated, with the tool *om_cgal_mesh_function* using the implicit function hemisphere called using the option *'-hr'* which is the hemisphere radius::

    ./tools/om_cgal_mesh_function -hr 0.87 -fs 0.5 -fd 0.05 -o northhemisphere.vtk

this generates the northern hemisphere, which we create a mirror image to create the southern hemisphere with matching vertices at their interface.

For more help on the cgal tools *om_cgal_mesh_function* see::

   ./tools/om_cgal_mesh_function -h
   ./tools/om_cgal_remesh -h
   ./tools/om_cgal_mesh_3Dlabeled_image -h
   ./tools/om_cgal_mesh_3Dlevelset_image -h

Using a tool such as Paraview_  we substract from these meshes their common interface (called *cut.vtk*), and merge all meshes into a single vtp file while naming these meshes::
 
  ./tools/om_mesh_to_vtp -i1 north.vtk -i2 south.vtk -i3 skull.vtk -i4 scalp.vtk -i5 cut.vtk -n1 "north" -n2 "south" -n3 "skull" -n4 "scalp" -n5 "cut" -o HeadNNc1.vtp
   

The files generated can easily be viewed using Paraview.

.. _CGAL: http://www.cgal.org/
.. _Paraview: http://www.paraview.org/
