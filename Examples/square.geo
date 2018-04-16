// Gmsh project created on Sun Mar 25 22:29:18 2018
SetFactory("OpenCASCADE");
//+
Point(1) = {-5, -5, 0, 1.0};
//+
Point(2) = {5, -5, 0, 1.0};
//+
Point(3) = {5, 5, 0, 1.0};
//+
Point(4) = {-5, 5, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Line Loop(1) = {1,2,3,4};
//+
Plane Surface(1) = {1};
