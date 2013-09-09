==========
= Models =
==========

Nested models:
^^^^^^^^^^^^^^
Head1, Head2, Head3 are 3 nested spheres as head models, with ratio 0.87, 0.92 and 1.

Head1 model contains 42 vertices per surface
Head2 model contains 162 vertices per surface
Head3 model contains 642 vertices per surface


Non-Nested models:
^^^^^^^^^^^^^^^^^^
- HeadNNa1, HeadNNa2, HeadNNa3, 
- HeadNNb1, HeadNNb2, HeadNNb3,
- HeadNNc1, HeadNNc2, HeadNNc3 are Non-Nested models.

- HeadNNaX is HeadX, where the inner sphere is splitted into two hemispheres: north & south. ( with X={1,2,3})
- HeadNNcX are the same models as HeadNNaX with a different(finer) discretization around the two hemispheres: north & south. ( with X={1,2,3})
         HeadNNc1 and HeadNNc3 were generated using the CGAL tool provided, while model HeadNNc2 is one of the used model in article of Jan Kybic. "Beyond nested volume...". ( with X={1,2,3})

- HeadNNbX is HeadX, where in the inner sphere were added two little spheres of radius 0.3 with center at (0., +/-0.4, 0.): spherenorth & spheresouth. ( with X={1,2,3})


================================
= Writting geom and cond files =
================================

\*.geom files:
^^^^^^^^^^^^^^

OpenMEEG definitions:
---------------------
"Mesh": a mesh is a collections of vertices and triangles all connected.
"Interface": an interface is a closed mesh.
"Half-space": an interface then defines a Inside and Outside. An Half-space is thus a space INSIDE or OUTSIDE an interface.
"Domain": a domain is a intersection of half-spaces.


namings:
--------

For ease of explanation, we consider the naming as follow:
- meshes: are named in lower case
- interfaces: in CamelCase
- domains: in UPPERCASE

but this convention is up to you.


A *.geom file* contains:

1. A header section: 
--------------------

|# Domain Description 1.1                             ||

which shows a version >= 1.1 since non-nested geometries

2. A (optional) MeshFile section: 
---------------------------------

|MeshFile Head.vtp                                    ||

It tells which file contains the geometry.
All meshes are in a single VTK/vtp file, where all polygons (triangles) have a string data associated
which tells the name of the surface it belongs. (these files can easily be opened with Paraview www.paraview.org, select some triangles-> Cell Label-> check Visible, see HeadNNa1/HeadNNa1.png for example )

3. An Interface section:
------------------------

Starting with the keyword "Interfaces", it follows the number of interfaces, and the keyword "Interface" or "Mesh".

In case no MeshFile section was found, "Mesh" must be the keyword. Else use "Interface".

|Interfaces 5 Interface                     ||        |Interfaces 3 Mesh
|                                           ||        |
|Interface North: +north +cut               ||        |Interface: "cortex.1.tri"
|Interface South: +south -cut               ||        |Interface: "skull.1.tri"
|Interface Cortex: north south              ||        |Interface: "scalp.1.tri"
|Interface Skull: skull                     ||        |
|Interface Scalp: scalp                     ||        |

Interface can be named or not (in both case). If not named, it will be automatically named '1', '2', ... following the order.

When using "Interface" keyword (on the left), meshes are orientable with a plus or minus sign, the overall interface (composed by several meshes or only one) should be consistently oriented (OpenMEEG will complain in other cases). 

3. A Domain section:
--------------------

|Domains 5                                  ||
|                                           ||
|Domain NORTH: -North                       ||
|Domain SOUTH: -South                       ||
|Domain SKULL: -Skull +Cortex               ||
|Domain SCALP: -Scalp +Skull                ||
|Domain AIR: +Scalp                         ||

Last, the definition of the domains using the previously defined interfaces.
A Domain is defined as being OUTSIDE or INSIDE certains interfaces.



\*.cond files:
^^^^^^^^^^^^^^

A *.cond file* contains:

1. A header section: 
--------------------

|# Properties Description 1.0 (Conductivities) ||

2. conductivities: 
------------------

|AIR         0.0                             ||
|NORTH       1                               ||
|SOUTH       1                               ||
|SKULL       0.0125                          ||
|SCALP       1                               ||
                                                
Just the domain names with there normalized conductivities.
