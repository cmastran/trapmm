$fn=32;
number_of_ribs= 16.000000;
plate_radius= 15.000000;
plate_thickness= 1.000000;
min_rib_width= 0.600000;
max_rib_width= 4.000000;
rib_radius= 6.708204;
rib_thickness= 0.200000;
rib_length= 8.775000;

module pbase() { 
	 polygon(points=[[-rib_length/2.0,min_rib_width/2.0],[rib_length/2,max_rib_width/2.0],[rib_length/2.0,-max_rib_width/2.0],[-rib_length/2.0,-min_rib_width/2.0]],paths=[[0,1,2,3]]); 
 }


module trapezoidal_rib() { 
	  linear_extrude(height=rib_thickness,center=true) { 
		 pbase(); 	 
 } 
 }
module trapezoidal_ribs() { 
 	 for (i = [0:number_of_ribs-1] ) { 
 		 angle = i*360.0/number_of_ribs; 
 		 rotate([0,0,angle]) { 
 			  translate([rib_radius,0,rib_thickness/2.0-0.005+plate_thickness/2.0]) trapezoidal_rib(); 
 		 } 
 	 } 
 } 
 

union() { 
 	 trapezoidal_ribs(); 
 	 cylinder(height=plate_thickness,r=plate_radius,center=true); 
 }
