//exact analytic solution


  double vertical_correction_min, vertical_correction_max, vertical_correction_c;   
  double r_ave_correction, C1, C2; 
  double max_edge_x_l, max_edge_x_u, max_edge_y_u, max_edge_y_l; 
  double min_edge_x_l, min_edge_x_u, min_edge_y_u, min_edge_y_l; 
  double xlim_l, xlim_u; //integral limits for partial intersections
  double xlim_mid; //middle integral limit for when cell fully contains beam

#ifdef TRASH
  for (k=0; k<ke; k++){
    for (j=js; j<je; j++){
      for (i=is; i<ie; i++){ //if the center position is in the path?
	//Centered angle of this solid angle bin
	if (k >= source_k_start && k <= source_k_end){ //only computing active solid angle bins
	  if (i == source_i && j == source_j)
	    realI[index] = 1.0; 
	  else if (i >= source_i && j >= source_j){ // exploiting the fact that all rays are in the first quadrant

	    if (xa1_b[k] == M_PI/2.0 || xa1_b[k] == 3*M_PI/2.0){ //vertical beam edge
	      vertical_correction_min = 1000000; // ensure that the source point (cone origin) is below the bottom corner and above top corner
	      slope_min = 0;//undefined slope 
	    }
	    else {
	      slope_min = tan(xa1_b[k]);
	      vertical_correction_min = 0.0; 
	    }
	    if (xa1_b[k+1] == M_PI/2.0 || xa1_b[k+1] == 3*M_PI/2.0){ //vertical beam edge
	      slope_max = 0; 
	      vertical_correction_max = 1000000;}
	    else{
	      slope_max = tan(xa1_b[k+1]);
	      vertical_correction_max = 0; 
	    }
	    if (xa1_b[k] == M_PI/2.0 || xa1[k] == 3*M_PI/2.0){ //vertical beam edge
	      slope_c = 0;  
	      vertical_correction_c = 1000000;}
	    else {
	      slope_c = tan(xa1[k]);
	      vertical_correction_c = 0;}	      
	  
	    dist = sqrt((x1[i] - source_x)*(x1[i] - source_x) + (x2[j] - source_y)*(x2[j] - source_y)); 
	    //Top beam edge linear plot locations
	    max_edge_x_l = source_y +vertical_correction_max + slope_max*(x1_b[i] - source_x); // y height of beam at x_{i-1/2}
	    max_edge_x_u = source_y +vertical_correction_max + slope_max*(x1_b[i+1] - source_x); // y height of beam at x_{i+1/2}
	    if (slope_max != 0){ //do not divide by zero 
	      max_edge_y_l = source_x +(x2_b[j] - source_y)/slope_max; // x position of beam edge when crossing bounday y_{j-1/2}
	      max_edge_y_u = source_x +(x2_b[j+1] - source_y)/slope_max; // x position of beam edge when crossing bounday y_{j-1/2}
	    }
	    else { // horizontal or vertical beam
	      if (xa1_b[k] == M_PI/2.0 || xa1[k] == 3*M_PI/2.0){//vertical
		max_edge_y_l = source_x;
		max_edge_y_u = source_x;  
	      }
	      else if (xa1_b[k] == 0){//horizontal beam is below all cells (never crosses any cell y boundaries 
		max_edge_y_l = 1000000;
		max_edge_y_u = 1000000;  
	      }
	    }
	    /* Compute y-heights of top left/bottom right corners of cell in question and corresponding beam boundary heights */
	    /* This trick only works if the beam is wholly in the right half of the plane */

	    //might be sufficient to simply check if one of the boundary rays lies inbetween y-edges or if cell is completely contained by beam: i.e. check positions of beam relative to both corners instead of vice versa
	    //debugging everything here:
	    if (i==source_i+1 && j==source_j){
	      printf("k = %d, cell within top beam = %d cell within bottom beam = %d\n",k,(x2_b[j+1] <= source_y +vertical_correction_max + slope_max*(x1_b[i] - source_x)),(x2_b[j] >= source_y -vertical_correction_min + slope_min*(x1_b[i+1] - source_x)));
	      printf("bottom edge height =%lf bottom, right corner beam height = %lf\n",x2_b[j],source_y -vertical_correction_min + slope_min*(x1_b[i+1] - source_x));
	      printf("top edge height =%lf top, left corner beam height = %lf\n",x2_b[j+1],source_y +vertical_correction_max + slope_max*(x1_b[i] - source_x)); }
	    //end debug

	    /* Cell fully contained within the beam path: both corners inside boundaries */ 
	    if ((x2_b[j+1] <= max_edge_x_l)// only need to check extreme corners due to first quadrant monotonicity
		&& (x2_b[j] >= min_edge_x_u)){
	      realI[index] = (dx1/2)*1.0/dist;
	      //try computing exact average. See notes for integral evaluation
	      C1 = (x1_b[i] - source_x)*(x1_b[i] - source_x); 
	      C2 = (x1_b[i+1] - source_x)*(x1_b[i+1] - source_x); 
	      r_ave_correction = 1.0/(4*dx2)*((x2_b[j+1] - source_y)*log((C2+pow((x2_b[j+1] - source_y),2))/(C1+pow((x2_b[j+1] - source_y),2))) + 2*sqrt(C2)*atan((x2_b[j+1] - source_y)/sqrt(C2)) - 2*sqrt(C1)*atan((x2_b[j+1] - source_y)/sqrt(C1))); 
	      printf("correction y_{j+1/2} term = %lf\n",r_ave_correction);
	      printf("C_{i-1/2} = %lf C_{i+1/2} = %lf\n",C1,C2); 
	      printf("atan(%lf / %lf ) = %lf\n",(x2_b[j+1] - source_y),sqrt(C2),atan((x2_b[j+1] - source_y)/sqrt(C2)));
	      r_ave_correction -= 1.0/(4*dx2)*((x2_b[j] - source_y)*log((C2+pow((x2_b[j] - source_y),2))/(C1+pow((x2_b[j] - source_y),2))) + 2*sqrt(C2)*atan((x2_b[j] - source_y)/sqrt(C2)) - 2*sqrt(C1)*atan((x2_b[j] - source_y)/sqrt(C1))); // should I use atan2?
	      printf("total correction term = %lf\n",r_ave_correction);
	      realI[index] += r_ave_correction; 
		printf("(%d,%d,%d) Contained cell average = %lf\n",i,j,k,realI[index]); 
	      if (i==source_i+1 && j==source_j){
		printf("k = %d fully contained\n",k);
	      }
	    }
	  
	    /* Top corner outside beam, bottom corner inside beam */
	    else if ((x2_b[j+1] >= max_edge_x_l)
		     && (x2_b[j] >= min_edge_x_u)
		     && (x2_b[j] <= max_edge_x_u)){
	      //first approach/approximation: if cell intersects beam at all, take the cell center distance from source to set the cell average intensity
	      //	      realI[index] = (dx1/2)*1.0/dist;
	      /* Second, exact analytic approach */
	      //Find where the top beam enters in x
	      if (max_edge_y_l > x1_b[i]){ //lower integral limit is not x_{i-1/2}
		//KYLE change this loop to ternary operator
		xlim_l = max_edge_y_l; 
	      }
	      else{
		xlim_l = x1_b[i]; 
	      }
	      //Find where the top beam exits cell in x
	      r_ave_correction = 0.0; 
	      if (max_edge_y_u < x1_b[i+1]){ //exits top y_{j+1/2} boundary
		// do two part integral
		xlim_u = max_edge_y_u; 
		// add on integral of remaining square region
		C1 = (xlim_u - source_x)*(xlim_u - source_x); 
		C2 = (x1_b[i+1] - source_x)*(x1_b[i+1] - source_x); 
		r_ave_correction = 1.0*(x1_b[i+1] - xlim_u)/2.0; 
		r_ave_correction += 1.0/(4*dx2)*((x2_b[j+1] - source_y)*log((C2+pow((x2_b[j+1] - source_y),2))/(C1+pow((x2_b[j+1] - source_y),2))) + 2*sqrt(C2)*atan((x2_b[j+1] - source_y)/sqrt(C2)) - 2*sqrt(C1)*atan((x2_b[j+1] - source_y)/sqrt(C1))); 
		r_ave_correction -= 1.0/(4*dx2)*((x2_b[j] - source_y)*log((C2+pow((x2_b[j] - source_y),2))/(C1+pow((x2_b[j] - source_y),2))) + 2*sqrt(C2)*atan((x2_b[j] - source_y)/sqrt(C2)) - 2*sqrt(C1)*atan((x2_b[j] - source_y)/sqrt(C1)));  
	      }
	      else{
		//do single integral with limits that depend on x
		xlim_u = x1_b[i+1]; 
	      }
	      //triangular region
	      //ensure referneces to slope_max never include the vertical lines
	      r_ave_correction += 0.5/dx2*((source_y - x2_b[j] - slope_max*source_x)*(xlim_u-xlim_l) +0.5*slope_max*(xlim_u*xlim_u - xlim_l*xlim_l)); 
	      C1 = (x2_b[j] - source_y)*(x2_b[j] - source_y);
	      r_ave_correction += 0.25/dx2*((xlim_u - source_x)*log(((1+slope_max*slope_max)*pow(xlim_u-source_x,2))/(C1 + pow(xlim_u-source_x,2))) - 2*sqrt(C1)*atan((xlim_u - source_x)/sqrt(C1)));
	      r_ave_correction -= 0.25/dx2*((xlim_l - source_x)*log(((1+slope_max*slope_max)*pow(xlim_l-source_x,2))/(C1 + pow(xlim_l-source_x,2))) - 2*sqrt(C1)*atan((xlim_l - source_x)/sqrt(C1)));
	      realI[index] = r_ave_correction; 
	    }
	    /* Top corner inside beam, bottom corner outside beam */
	    else if ((x2_b[j+1] <= max_edge_x_l)
		     && (x2_b[j] <= min_edge_x_u)
		     && (x2_b[j+1] >= min_edge_x_l)){
	      //first approach/approximation: if cell intersects beam at all, take the cell center distance from source to set the cell average intensity
	      //       realI[index] = (dx1/2)*1.0/dist; 	      
	      /* Second, exact analytic approach */
	      //Find where the bottom beam exits cell in x
	      if (min_edge_y_u < x1_b[i+1]){ //exits top y_{j+1/2} boundary
		xlim_u = min_edge_y_u; 
	      }
	      else{
		//do single integral with limits that depend on x
		xlim_u = x1_b[i+1]; 
	      }
	      //Find where the bottom beam edge enters cell in x
	      r_ave_correction = 0.0; 
	      if (min_edge_y_l > x1_b[i]){ //lower integral limit is not x_{i-1/2}
		//KYLE change this loop to ternary operator
		// do two part integral
		xlim_l = min_edge_y_l; 
		// add on integral of remaining square region
		C1 = (x1_b[i]  - source_x)*(x1_b[i] - source_x);//C1,C2 contain all references to x limits 
		C2 = (xlim_l - source_x)*(xlim_l - source_x); 
		r_ave_correction = 1.0*(xlim_l - x1_b[i])/2.0; 
		r_ave_correction += 1.0/(4*dx2)*((x2_b[j+1] - source_y)*log((C2+pow((x2_b[j+1] - source_y),2))/(C1+pow((x2_b[j+1] - source_y),2))) + 2*sqrt(C2)*atan((x2_b[j+1] - source_y)/sqrt(C2)) - 2*sqrt(C1)*atan((x2_b[j+1] - source_y)/sqrt(C1))); 
		r_ave_correction -= 1.0/(4*dx2)*((x2_b[j] - source_y)*log((C2+pow((x2_b[j] - source_y),2))/(C1+pow((x2_b[j] - source_y),2))) + 2*sqrt(C2)*atan((x2_b[j] - source_y)/sqrt(C2)) - 2*sqrt(C1)*atan((x2_b[j] - source_y)/sqrt(C1)));  
	      }
	      else{
		//do one integral
		xlim_l = x1_b[i]; 
	      }
	      //triangular region
	      //ensure referneces to slope_min never include the vertical lines
	      r_ave_correction += 0.5/dx2*((x2_b[j+1] - source_y  + slope_min*source_x)*(xlim_u-xlim_l) - 0.5*slope_min*(xlim_u*xlim_u - xlim_l*xlim_l)); 
	      C2 = (x2_b[j+1] - source_y)*(x2_b[j+1] - source_y);
	      r_ave_correction += 0.25/dx2*((xlim_u - source_x)*log((C2 + pow(xlim_u-source_x,2))/((1+slope_min*slope_min)*pow(xlim_u-source_x,2))) - 2*sqrt(C2)*atan((xlim_u - source_x)/sqrt(C2)));
	      r_ave_correction -= 0.25/dx2*((xlim_l - source_x)*log((C2 + pow(xlim_l-source_x,2))/((1+slope_min*slope_min)*pow(xlim_l-source_x,2))) - 2*sqrt(C2)*atan((xlim_l - source_x)/sqrt(C2)));
	      realI[index] = r_ave_correction; 
	    }
	      
	    /* Cell fully contains beam */
	    else if ((x2_b[j+1] >= max_edge_x_l)
		     && (x2_b[j+1] >= min_edge_x_l)
		     && (x2_b[j] <= max_edge_x_u)
		     && (x2_b[j] <= max_edge_x_l)){
	      //realI[index] = (dx1/2)*1.0/dist;
	      //top edge of beam must enter and exit first (in x positions at y_j-1/2 and +1/2 respectively)

	      //Exits:
	      //top beam exits through top boundary
	      if (max_edge_y_u < x1_b[i+1]){
		//Exit case #1: both beams exit through top y_{j+1/2} boundary
		if (min_edge_y_u < x1_b[i+1]){
		  xlim_u = min_edge_y_u; 
		  //Beam entrances:
		  //Entrance a) top edge of beam enters through lower y_{j-1/2} boundary, so bottom edge must also enter through here
		  if (max_edge_y_l > x1_b[i]){
		    xlim_l = max_edge_y_l; 
		    xlim_mid = min_edge_y_l;
		  }		
		  else{
		    //top edge of beam enters through left x_{i-1/2} boundary. Bottom edge can enter here or bottom boundary
		    xlim_l = x1_b[i];
		    //Entrance b) bottom edge enters through bottom boundary 
		    if (min_edge_y_l > x1_b[i]){ //two integrals
		      xlim_mid = min_edge_y_l;
		    }
		    else{ //Entrance c) both enter through left boundary
		      // two integrals
		      xlim_mid = max_edge_y_u; 
		    }		    
		  }		  
		}
		else{ //Exit case #2: top beam exits through top boundary, bottom beam exits through right boundary
		  xlim_u = x1_b[i+1];
		  //Beam entrances:
		  //Entrance a) top edge of beam enters through lower y_{j-1/2} boundary, so bottom edge must also enter through here
		  if (max_edge_y_l > x1_b[i]){ //double check this case KYLE
		    xlim_l = max_edge_y_l; 
		    xlim_mid = min_edge_y_l; //can these be swapped??
		    xlim_mid2 = max_edge_y_u; 
		  }		
		  else{
		    //top edge of beam enters through left x_{i-1/2} boundary. Bottom edge can enter here or bottom boundary
		    xlim_l = x1_b[i];
		    //Entrance b) bottom edge enters through bottom boundary 
		    if (min_edge_y_l > x1_b[i]){ 
		      //3 integrals
		      xlim_mid = max_edge_y_u; 
		      xlim_mid2 = min_edge_y_l;
		    }
		    else{ //Entrance c) both enter through left boundary
		      // two integrals
		      xlim_mid = max_edge_y_u; 
		    }		    
		  }
	      }
	      else{
		//Exit case #3: both beams exit through x_{i+1/2} boundary
		xlim_u = x1_b[i+1]; 
		  //Beam entrances:
		  //Entrance a) top edge of beam enters through lower y_{j-1/2} boundary, so bottom edge must also enter through here
		  if (max_edge_y_l > x1_b[i]){
		    xlim_l = max_edge_y_l; 
		    xlim_mid = min_edge_y_l;
		  }		
		  else{
		    //top edge of beam enters through left x_{i-1/2} boundary. Bottom edge can enter here or bottom boundary
		    xlim_l = x1_b[i];
		    //Entrance b) bottom edge enters through bottom boundary 
		    if (min_edge_y_l > x1_b[i]){ //two integrals
		      xlim_mid = min_edge_y_l;
		    }
		    else{ //Entrance c) both enter through left boundary
		      // two integrals
		      xlim_mid = max_edge_y_u; 
		    }		    
		  }		  
	      }
	      }
	    }
	  }
	      else{
	      realI[index] = 0.0; 
	      }	    
	    }
	  }
	}
	else{
	  realI[index] = 0.0; 
	}
	index++; 
      }
    }
  }
#endif //TRASH
