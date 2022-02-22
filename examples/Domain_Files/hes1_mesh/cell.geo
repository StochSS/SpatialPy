
// Gmsh project created on Wed Jun 11 15:38:08 2014

lc = 0.9; //characteristic length
cr = 7.5; //cell radius
nr = 0.4; //nucleus size ratio

Point(1) = {0, 0, 0, lc};

//Cell
Point(2) = {cr, 0, 0, lc};
Point(3) = {0,cr, 0, lc};
Point(4) = {-cr, 0, 0, lc};
Point(5) = {0,-cr, 0, lc};
Circle(1) = {2, 1, 3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};
Rotate {{1,0,0},{0,0,0},Pi/2} {
  Duplicata{Line{1,2,3,4};}
}
Line Loop(9) = {8, 1, 2, 7};
Ruled Surface(10) = {9} In Sphere {1};
Line Loop(11) = {2, -6, -5, 1};
Ruled Surface(12) = {11} In Sphere {1};
Line Loop(13) = {6, 3, 4, 5};
Ruled Surface(14) = {13} In Sphere {1};
Line Loop(15) = {7, 8, -4, -3};
Ruled Surface(16) = {15} In Sphere {1};
Surface Loop(17) = {12, 10, 16, 14};
Volume(18) = {17};
//Nucleus
Dilate { {0,0,0}, nr } {
  Duplicata{Volume{18};}
}
Delete { Volume{18};}

Surface Loop(36) = {25, 30, 35, 20};
Volume(37) = {17, 36};
//Physical Point(0) = {1};
Physical Surface(0) = {20, 35, 30, 25};
Physical Surface(1) = {10, 12, 14, 16};
Physical Volume(1) = {19};
Physical Volume(2) = {37};
