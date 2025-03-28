// Nicholas Rathmann, 2024

scale=1;
scale=1.5;

lc=8.0e3 * scale;
lc2=2*lc;

x0 = -900e3;
x1 = +600e3;
y0 = -1400e3;
y1 = -400e3;


Point(1)={-550e3,-1225e3, 0.0,lc}; //1
Point(2)={-650e3,-1100e3, 0.0,lc};
Point(3)={x0,-1000e3, 0.0,lc};
Point(4)={x0,y1, 0.0,lc};
Point(5)={x1,y1, 0.0,lc};
Point(6)={x1,-1250e3, 0.0,lc};
// free bnd
Point(7)={300e3,-1300e3, 0.0,lc};
Point(8)={200e3,-1300e3, 0.0,lc};
Point(9)={-350e3,-1150e3, 0.0,lc};

Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,5};
Line(5)={5,6};
Line(6)={6,7};
Line(7)={7,8};
Line(8)={8,9};
Line(9)={9,1};

Line Loop(5)={1,2,3,4,5,6,7,8,9};

Plane Surface(10)={5};

Physical Surface(11)={10};

Physical Line (1) = {1,2,3,4,5,6}; // isotropic bnd
Physical Line (2) = {7,8,9}; // free bnd

// --------------------

Point(1001)={00e3,-0.4e6,0,lc};
Point(1002)={600e3,-0.6e6,0,lc};
Line(50) = {1001,1002};

Point(1003)={600e3,-0.4e6,0,lc};
Point(1004)={600e3,-1.3e6,0,lc};
Line(51) = {1003,1004};

Point(1005)={-500e3,-0.65e6,0,lc};
Point(1006)={-350e3,-0.75e6,0,lc};
Line(52) = {1005,1006};

Point(1007)={0e3,-1e6,0,lc};
Point(1008)={0e3,-1.3e6,0,lc};
Line(53) = {1007,1008};

Point(1009)={x0-40e3,y0,0,lc};
Point(1010)={x0-40e3,y1,0,lc};
Line(54) = {1009,1010};

// Say we would like to obtain mesh elements with size lc/30 near curve 2 and
// point 5, and size lc elsewhere. To achieve this, we can use two fields:
// "Distance", and "Threshold". We first define a Distance field (`Field[1]') on
// points 5 and on curve 2. This field returns the distance to point 5 and to
// (100 equidistant points on) curve 2.
Field[1] = Distance;
Field[1].CurvesList = {50,51,52,53,54};

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
Field[2].DistMin = 00e3;
Field[2].DistMax = 150e3;

Background Field = 2;

// Don't extend the elements sizes from the boundary inside the domain
Mesh.CharacteristicLengthExtendFromBoundary = 0;
