// This variable is a characteristic length and it controles the mesh size around a Point.

// It is possible to specify more than one variable for this purpose.
// 1.5 1.25 1 0.75 0.4 0.2 0.1 
	cl1 = 0.1*1.5;




// Points contains the x, y and z coordinate and the characteristic length of the Point.


	Point(1) = {0,0,0,cl1};

	Point(2) = {0,0.5,0,cl1};

	Point(3) = {-0.5,0,0,cl1};



// A Line is basically a connection between two Points. A good practice is to connect the

// Points in a counter-clockwise fashion.


	Line(1) = {1,2};
	Circle(3) = {2,1,3};
	
Line(2) = {3,1};
	


// A Line Loop is the connection of Lines that defines an area. Again it is good practice

// to do this in a counter-clockwise fashion.


	Line Loop(1) = {1,3,2};



// From the Line Loop it is now possible to create a surface, in this case a Plane Surface.


	Plane Surface(1) = {1};



// From the Plane Surface a Physical Surface is generated, this makes is possible to only

// save elements which are defined on the area specified by the Line Loop.


	Physical Surface(1) = {1};