// *** Mesh of Ross, Antarctica, for steady state SSA CPO solver ***
// Lines 3 to 6 below *must* be the bounding box coordinates of the model domain (to be read by python code)
x0=-900e3;
y0=-1320e3;
x1=+600e3;
y1=-400e3;

// Mesh resolution
//scale=1;
scale=1.5;
lc=8.0e3 * scale;
lc2=2*lc;

// Mesh corners
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

// Connect corners
Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,5};
Line(5)={5,6};
Line(6)={6,7};
Line(7)={7,8};
Line(8)={8,9};
Line(9)={9,1};
Line Loop(50)={1,2,3,4,5,6,7,8,9};

// Physical lines, etc.
Physical Line(1)={1,2,3,4,5,6}; // isotropic boundary on grounded ice
Physical Line(2)={7,8,9};       // free boundary on shelf edge
Plane Surface(10)={50};
Physical Surface(11)={10};

// Refine mesh
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

Field[1] = Distance;
Field[1].CurvesList = {50,51,52,53,54};
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lc2;
Field[2].SizeMax = lc;
Field[2].DistMin = 00e3;
Field[2].DistMax = 150e3;
Background Field = 2;
Mesh.CharacteristicLengthExtendFromBoundary = 0;
