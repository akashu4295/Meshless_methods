SetFactory("OpenCASCADE");

Circle(1) = {0, 0, 0, 1, 0, 2*Pi};
Circle(2) = {0, 0, 0, 2, 0, 2*Pi};

Curve Loop(1) = {2};
Curve Loop(2) = {1};

Plane Surface(1) = {1, 2};

Physical Curve("innerwall", 3) = {1};
Physical Curve("outerwall", 4) = {2};
Physical Surface("fluid", 5) = {1};

Mesh.MeshSizeMax = 0.04;
