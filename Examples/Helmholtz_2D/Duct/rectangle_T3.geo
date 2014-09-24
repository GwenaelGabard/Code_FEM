H =    2.00000;
L =    2.00000;
MeshLength =    0.05000;




 //////////////////////////////////////////////////////////////////
 
Mesh.ElementOrder = 1;

// Physical surface domain
DOMAIN_1   = 1000;

BOUNDARY_1001 = 1001;  
BOUNDARY_1002 = 1002; 
BOUNDARY_1003 = 1003; 
BOUNDARY_1004 = 1004; 


Point(4) = {0, 0, 0, MeshLength};
Point(1) = {L, 0, 0, MeshLength};
Point(2) = {L, H, 0, MeshLength};
Point(3) = {0, H, 0, MeshLength};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(9) = {4,1,2,3};
Plane Surface(10) = {9};

// Define the physical groups to apply the boundary conditions
Physical Line(BOUNDARY_1001) = {4};
Physical Line(BOUNDARY_1002) = {1};
Physical Line(BOUNDARY_1003) = {2};
Physical Line(BOUNDARY_1004) = {3};

// Physical domain
Physical Surface(DOMAIN_1) = {10};
