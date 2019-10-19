// This variable is a characteristic length and it controles the mesh size around a Point.
// It is possible to specify more than one variable for this purpose.
cl1 = 0.02288/10;
cl2 = 0.02288/10;
// Points contains the x, y and z coordinate and the characteristic length of the Point.
Point(1) = {0,0,0,cl1};
Point(2) = {0.02288,0,0,cl1};
Point(3) = {0.02288,0.02288,0,cl2};
Point(4) = {0,0.02288,0,cl2};

// A Line is basically a connection between two Points. A good practice is to connect the
// Points in a counter-clockwise fashion.
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

// A Line Loop is the connection of Lines that defines an area. Again it is good practice
// to do this in a counter-clockwise fashion.
Line Loop(1) = {1,2,3,4};

// From the Line Loop it is now possible to create a surface, in this case a Plane Surface.
Plane Surface(1) = {1};

Physical Surface("rve0") = {1};
