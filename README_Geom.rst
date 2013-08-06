OpenMEEG: Taking into account non-nested geometries
===================================================

Definitions:
------------
-Vertex: 3D vector of real                                                     Vect3
 vertices: collection of vertex                                                Vect3* vertices   OR    Vertex vertices
 
-Triangle: vector of 3 indices                                                 Triangle

-Mesh: collection of triangles                                                 Triangle* triangles

-Interface: A closed-connexed shape, e.g a mesh homeomorphe to a sphere.       vector<Mesh&> interface

-OrientedInterface: An Oriented (with normals) Interface                       vector<pair< interface, vector<bool> > > Ointerface
                    the boolean tell wheater the direct product is well 
                    oriented or not with regard to the interior domain.

-Domain: 3D volume separated by 2-(or more) interfaces (or 1 OrientedInterface, e.g the Air domain)

(words we do not used: layers, face, cells, surface )

Features:
---------
We would like to get a description of the geometry such that this API gives us access to the following features::
    - Given a domain ID::
            - define the interfaces enclosing this domain: RETURN std::vector<Interface> 
            - define the domain conductivity, and if a reference to a triangle is also given, gives the conductivity over the triangle

    - Given a 3D point, gives the ID of the domain which contains the point

    - Given an Interface, able to iterate over the triangles/points of this interface

    - Acquire geometry description through .geom and .cond files as well as more general VTK format.
    - Read/write mesh files in different format
    - Test wheather or not a mesh is closed, self intersecting,...
    - Determine gradient of P1 functions, and surface of triangles
    - ...

    - Creates Interface from meshes


Examples of Non-Nested Geometry:
--------------------------------
     ____________              _____________
    /            \            /             \
   /              \          /               \
   |   ______     |          |   ___   ___   |
   |  /   |  \    |          |  /   |  |  \  |
   |  | 1 | 2 |   |    OR    |  | 1 |  |2  | |
   |  \___|__/    |          |  \___|  |__/  |
   |              |          |               |
   \      3      /           \      3       /
    \___________/             \____________/
                          

VTK based format::
==================                                       _
                                                ||      /    Vertices:                             ||
                                                ||      |        x1 y1 z1                          ||
                                                ||      |        x2 y2 z2                          ||
                                                ||      |        .. .. ..                          ||
                                                ||      |        xN yN zN                          ||
                                                ||      |    Triangles:                            ||
                                                ||      |        p1 p2 p3                          ||
                                                ||      |        p3 p7 p6                          ||
                                                ||      |        p7 p8 p9                          ||
                                                ||      |        .. .. ..                          ||
                                                ||      |    Meshes:                               ||
                                                ||      |       m0:                                ||
                                                ||      |         t0                               ||
                                                ||      |         t1                               ||
                                                ||  VTK |         ..                               ||
                                                ||      |       m1:                                ||
                                                ||      |         t12                              ||
                                                ||      |         t37                              ||
                                                ||      |         t78                              ||
                                                ||      |         ...                              ||
                                                ||      |       ..:                                ||
                                                ||      \_        ...                              ||
# Domain Description 1.0                        ||      /    # Domain Description 1.1              ||      /    # Domain Description 1.1
                                                ||      |                                          ||      |
Interfaces 3 Mesh                               ||      |    MeshFile file.vtp                     ||      |                 
                                                ||      |                                          ||      |    
"skull.1.tri"                                   ||      |    Interfaces 3 NamedInterface           ||      |                   
"cortex.1.tri"                                  ||      |       i0:  m0 m1 m2 ..                   ||      |                    
"scalp.1.tri"                                   ||      |       i1:  m1 m8 ..                      ||      |                    
                                                ||      |       ..:  .. ..                         ||      |            
Domains 4                                       ||      |    Domains 4                             ||      |    Interfaces 3 NamedMeshes   
                                                ||      |        Domain GrayMatter_Left:  -i0      ||      |       i0: skull.1.tri
Domain Scalp 1 -3                               ||      |        Domain GrayMatter_Right: -i1      ||      |       i1: cortex.1.tri
Domain Brain -2                                 ||      |        Domain CSF:        i0 i1 -i2      ||      |       i2: scalp.1.tri
Domain Air 3                                    ||      |        Domain Air:        i2             ||      |     
                                                || GEOM |                                          || GEOM |    Domains 4                          
Domain Skull 2 -1                               ||      |                                          ||      |        Domain Scalp -i1
                                                ||      |                                          ||      |        Domain Brain -i2
                                                ||      |                                          ||      |        Domain Air    i0 i1 -i2
                                                ||      |                                          ||      |        
                                                ||      |                                          ||      |        Domain Skull 2 -1 
                                                ||      |                                          ||      |    
                                                ||      |                                          ||      |    
                                                ||      |                                          ||      |        
                                                ||      \_                                         ||      \_     

Design of the reconstruction:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    Remove the current files mesh3, MeshReader, MeshDescription/*
    Add files from (Odyssee++/)FMesh, like Face.H, Triangle.H ... ??
    Replace the current geometry class by something for general


    class Reader
    class geometry, which construct domains
    TODO: specify that 2 meshes can have only one domain shared or zero, if they have 2 shared domains, then they are the same mesh.


Reading process of the geom file:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
read header: ()
    Domain Description: version 1.0
if (version == 1.0)
    read Interfaces nb_meshes Mesh
        for i in nb_meshes
            [p,t] = meshReader(mesh%i)
            points.pushback_ (p, from mesh %i)
            meshes.pushback_ (t, from mesh %i)
         end
         for i in meshes
             interfaces.pushback_ (&meshes[i])
         end
    
    read Domains nb_domains
        for i in nb_domains
            domains.pushback_ (name%i, +/- interface )
        end
else (version == 2.0) 
    read VTK file 
        points.pushback_ (VTK)
        meshes.pushback_ (VTK)
    read Interfaces nb_interfaces
        for i in nb_interfaces
            interfaces.pushback_ (&meshes[j])
        end
    read Domains nb_domains
        for i in nb_domains
            domains.pushback_ (name%i, +/- interface )
        end


TODO
^^^^

- Do we let the interfaces in the eometry class, or is it included in each Domain ?
- Do we keep a Reference on the Triangle class or * is ok ?  
- Remove 'inverse' in header



Questions Theo:

- Dans le geom:
* Mot-clef: Interface pour les Interfaces ? comme Domain/Domains
* Mettre un ":" apres le nom du domaine pour delimiter avec les noms des interfaces

- Dans geometry:
  * as t-on besoin de stocker les interfaces ? sachants qu'elles se trouvent en ce moment dans les domains.
  * Est-ce qu'une ionterface est une collection de paire de mesh, bool, pour dire si oui/non la mesh est orienté correctement ?
  * 

- VTK:
  * J'oblige à avoir VTK ? ou est-ce que ce sera que optionnel et avec des avantages en + (non nested) ?



  TODO:

  doxygen avec ///
  regarder la memoire et temps de calculs
  oprofile

  DONE:
  :tabdo %g/ &[a-zA-Z,]/s/ &/\& /gc
  (
  include <> instead <>
  ++i instead i++
  )
