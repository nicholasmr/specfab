// Nicholas Rathmann, 2024-

scale=1;
scale=1.7;

lc=2.5e3 * scale;
lc2=lc/2;

x0=-1710e3;
x1=-1425e3;
y1=-50e3;

Point(1)={x0,-304e3, 0.0,lc};
Point(2)={x0,y1, 0.0,lc};
Point(3)={x1,-50e3, 0.0,lc};
Point(4)={x1,-275e3, 0.0,lc};
Point(5)={-1520e3,-335e3, 0.0,lc};

Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,5};
Line(5)={5,1};

Line Loop(5)={1,2,3,4,5};

Physical Line (1) = {1,2,3,4}; // isotropic bnd
Physical Line (2) = {5}; // free bnd

Plane Surface(10)={5};

Physical Surface(11)={10};


// --------------------

Point(1001)={-1625e3,-350e3,0,lc};
Point(1002)={-1580e3,-230e3,0,lc};
Line(50) = {1001,1002};


// Say we would like to obtain mesh elements with size lc/30 near curve 2 and
// point 5, and size lc elsewhere. To achieve this, we can use two fields:
// "Distance", and "Threshold". We first define a Distance field (`Field[1]') on
// points 5 and on curve 2. This field returns the distance to point 5 and to
// (100 equidistant points on) curve 2.
Field[1] = Distance;
Field[1].CurvesList = {50};

// We then define a `Threshold' field, which uses the return value of the
// `Distance' field 1 in order to define a simple change in element size
// depending on the computed distances
//
// SizeMax -                     /------------------
//                              /
//                             /
//                            /
// SizeMin -o----------------/
//          |                |    |
//        Point         DistMin  DistMax

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lc2;
Field[2].SizeMax = lc;
Field[2].DistMin = 30e3;
Field[2].DistMax = 70e3;

Background Field = 2;

// Don't extend the elements sizes from the boundary inside the domain
Mesh.CharacteristicLengthExtendFromBoundary = 0;

