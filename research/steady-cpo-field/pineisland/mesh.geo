// Nicholas Rathmann, 2024-
// Lines 3 to 6 below *must* be the bounding-box coordinates of the model domain (to read by python code)
x0=-1710e3;
y0=-335e3;
x1=-1425e3;
y1=-50e3;

// Mesh resolution
lc=5e3;
lc2=lc/2;

// Mesh corners
Point(1)={x0,-304e3, 0.0,lc};
Point(2)={x0,y1, 0.0,lc};
Point(3)={x1,-50e3, 0.0,lc};
Point(4)={x1,-275e3, 0.0,lc};
Point(5)={-1520e3,-335e3, 0.0,lc};

// Connect corners
Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,5};
Line(5)={5,1};
Line Loop(50)={1,2,3,4,5};

// Physical lines, etc.
Physical Line (1)={1,2,3,4}; // isotropic boundary on grounded ice
Physical Line (2)={5};       // free boundary on shelf edge
Plane Surface(10)={50};
Physical Surface(11)={10};

// Refine the mesh in the neighborhood of a line that passes through flow trunk

Point(1001)={-1625e3,-350e3,0,lc};
Point(1002)={-1580e3,-230e3,0,lc};
Line(50)={1001,1002};

Field[1]=Distance;
Field[1].CurvesList={50};
Field[2]=Threshold;
Field[2].InField=1;
Field[2].SizeMin=lc2;
Field[2].SizeMax=lc;
Field[2].DistMin=30e3;
Field[2].DistMax=70e3;
Background Field=2;
Mesh.CharacteristicLengthExtendFromBoundary=0; // Don't extend the elements sizes from the boundary inside the domain
