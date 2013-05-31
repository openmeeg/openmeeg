OpenMEEG: Taking into account non-nested geometries
===================================================

Definitions:
------------
-Domain: Volume separated between interfaces

-Interface: A closed-connexed shape, e.g a mesh homeomorphe to a sphere.

-Mesh: collection of points and cells (triangles)


Features:
---------
We would like to get a description of the geometry such that this API gives us access to the following features::
    - Given a domain ID, define the interfaces enclosing this domain
    - Given a point, give the Id of the domain where this point lie
    - Able to iterate over the triangles/points of an interface
    - Read/write mesh files in different format
    - Acquire geometry description through .geom and .cond files as well as more general VTK format.
    - Test wheather or not a mesh is closed, self intersecting,...
    - Determine gradient of P1 functions, and surface of triangles
    - ...


Example of Non-Nested Geometry:
-------------------------------
     ____________
    /            \
   /              \
   |   ___ ___    |
   |  /   ||  \   |
   |  \_1_||_2/   |
   |              |
   \      3      /
    \___________/


VTK based format:
-----------------
Vertices::
    x1 y1 z1
    x2 y2 z2
    .. .. ..
    xN yN zN
Triangles::
    p1 p2 p3
    p3 p7 p6
    p7 p8 p9
    .. .. ..
    .. .. ..
Interfaces::
   i0: (GrayMatter_Left)
        t1
        t2
        t3
        t4
        ..
   i1: (GrayMatter_Right)
        t1
        t2
        t9
        t14
        ..
   i2: (Inner Skull)
        ..
        ..
   ..:
Domains(??? or in the .geom file)::
    GrayMatter_Left:    -i0
    GrayMatter_Right:   -i1
    CSF:                 i0 i1 -i2
    ..............................



Design of the reconstruction:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    Remove the current files mesh3, MeshReader, MeshDescription/*
    Add files from (Odyssee++/)FMesh, like Face.H, Triangle.H ... ??
    Replace the current geometry class by something for general

