TriangleMeshDistance
====================

Header only, single file, simple and efficient C++11 library to compute the signed distance function (SDF) to a triangle mesh.

The distance computation to the triangle collection is accelerated with a sphere bounding volume hierarchy. The sign of the distance is resolved with the method presented in *"Generating Signed Distance Fields From Triangle Meshes"* by Bærentzen, Andreas & Aanæs, Henrik. (2002).

Assuming triangle normals point outwards from the enclosed volume, the sign of the distance will be positive outside and negative inside.

https://github.com/InteractiveComputerGraphics/TriangleMeshDistance

License: MIT

Copyright (c) 2021 Jose Antonio Fernandez Fernandez
