#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "visit_writer.h"
#include "gsl/gsl_integration.h"

//code modifications to-do:
/* 1) Github version control
2) eliminate most preprocessor directives
3) Condense debug statements 
4) create variables to probe specific elements at specific timesteps (see 3)
5) is cfl condition 1/2????
6) reconcile error reporting (negative I, etc) with other codes
7) more status /progress print statements

8) mark all persistent errors (that may exist in other codes), with ERROR
9) fix 2D, 3D array allocation 
9) delete old functions
10) increase modularity. Add header file
 candidates for new functions:
- zero out realI, realJ
- copy I to realI, duplicate realJ
- analytic solution
-
11) really clean up flux loop
- rename variables to appropriate quantities: flux_limiter --> I_edge
12) compute numerical fluxes_rad
*/

/* Select solver options */
#define X1_PERIODIC 0
#define X2_PERIODIC 0
#define AUTO_TIMESTEP 1 //flag for automatically setting dt such that max{cfl_array} = CFL
#define CFL 0.8 //if set to 1.0, 1D constant advection along grid is exact. 
#define OUTPUT_INTERVAL 1 //choose timestep interval to dump simulation data to rte-n.vtk. 0 for only last step, 1 for every step
#undef SECOND_ORDER //flag to turn on van Leer flux limiting
#define ANALYTIC_SOLUTION //in addition to timestepping, output steady state solution in analytic-rte.vtk


#define V_EXACT //ifdef, then exact integral over boundary velocity is computed. Else, midpoint approximation is used
#undef DEBUG //define to probe a cell's relevant values at indices assigned below
#define TEST5
#undef SOLID_ANGLE //only works with first order. recompute the edge waves based on integral average of positive solid angle 
#undef FACE_INTERPOLATE //only works with #undef SECOND_ORDER, #undef SOLID_ANGLE

double X_physical(double, double);
double Y_physical(double, double);
double initial_condition(double x, double y, double xa1);
double bc_x1i(double x, double y, double xa1, double t);
double bc_x1f(double x, double y, double xa1, double t);
double bc_x2i(double x, double y, double xa1, double t);
double bc_x2f(double x, double y, double xa1, double t);
double flux_PLM(double ds,double *imu);
float find_max(float a[], int n);
float find_min(float a[], int n); 
float sum(float a[], int n);
double gaussian(double x_0, double y_0,double x,double y);
int uniform_angles2D(int N, double phi_z, double *pw, double **mu, double **mu_b, double *xa1,double *xa1_b);
double **allocate_2D_contiguous(double *a, int n1, int n2);
double ***allocate_3D_contiguous(double *a, int n1, int n2, int n3);
double find_max_double(double a[], int n);
float bc_x1f_polar(double phi, double xa1, double t);
double flux_PLM_athena(double r[3], int dir, double dt, double ds, double vel, double imu[3]);
double flux_PLM_athena_debug(double r[3], int dir, double dt, double ds, double vel, double imu[3]);
double bc_interior(double r, double phi, double xa1, double t);
void gaussianelim(double **A, double *b, double *x, int n, int pivot);
int permutation(int i, int j, int k, int ** pl, int np);
int bruls_angles2D(int N, double phi_z, double *pw, double **mu, double **mu_b, double *xa1,double *xa1_b);
int cmpfunc (const void * a, const void * b);

float analytic_solution(double x,double y, double r_max,double angle_c);

//functions for solid angle numerical quadrature
double max_sine (double x, void *params);
double min_sine (double x, void *params);
double max_cosine (double x, void *params);
double min_cosine (double x, void *params);

int main(int argc, char **argv){
  int i,j,k,l,n; 
  int index=0; 
  int indexJ=0;   /* for computing zeroth angular moment of the radiation field*/
  int next, next2, prev, prev2; //for indexing third/angular dimension since ghost cells are not used in that dimension
  int nsteps = 50;
  double dt = 0.02; //initial timestep. Auto adjusted #ifdef AUTO_TIMESTEP

  /*----------------------------------------*/
  /*         Initialize Spatial Mesh        */
  /*----------------------------------------*/
  printf("Initializing mesh.....\n");
  /* Computational (2D polar) grid coordinates */
  int nx1 = 128;
  int nx2 = 128; 

  /* Angular parameter */
  int xa1_uniform = 1; //boolean for switching between uniform discretization and Bruls discretization
  int N_bruls; //analogous to N in MATLAB source code. must be even and satisfy N_bruls <= 12
  int nxa1;
  double *dxa1;  //angular width of cell 
  if (xa1_uniform){
    nxa1 = 4; 
  } 
  else {
    N_bruls = 4;
    nxa1 = N_bruls*(N_bruls+2)/2;
  }
  dxa1 = malloc(sizeof(double)*nxa1); //dxa1[] is an array to allow for nonuniform option
  //polar angle relative to 2D spatial domain at which photons propagate in z. phi_z=0 ---> normal to surface
  double phi_z = M_PI/2; 

  //TEST5 parameters: these are the active angular bins (inclusive) that point source emits radiation
  int source_k_start = 0;
  int source_k_end = 3;

  /* DEBUG parameters: set these to select a cell and timestep at which to print out fluxes and limited slopes, etc */
#ifdef DEBUG
  int n_debug = 1;
  int k_debug = 0;
  int j_debug = 66;
  int i_debug = 67;
#endif 

  int num_ghost = 2; //number of ghost cells on both boundaries of each dimension
  /* num_ghost = 1 for piecewise constant reconstruction. 
     num_ghost = 2 for piecewise linear reconstruction 
     num_ghost = 3 for piecewise parabolic reconstruction */
  int nx1_r = nx1; //number of non-ghost cells in x
  int nx2_r = nx2; //number of non-ghost cells in y
  nx1 += 2*num_ghost; //total number of cells in x
  nx2 += 2*num_ghost; //total number of cells in y

  /* Non-ghost indices, r.f. ATHENA convention */
  // Use these to step through all physical loops
  int is = num_ghost;
  int ie= is+nx1_r; 
  int js = num_ghost;
  int je = js+nx2_r;
  int ks = 0;
  int ke = nxa1;

  /* CONVENTION: the above ranges refer to the non-ghost cells. However, the ghost cells may have real coordinate interpretations (i.e. they may represent actual cells right outside domain). In some circumstances this means that we must be sure that the number of ghost cells makes sense with the range of coordinates (more directed towards polar coordinates radius). */

  /* CONVENTION: all mesh structures (zonal and nodal) will be nx1 x nx2 x nxa1, although we may not fill the ghost entries with anything meaningful. This is to standardize the loop indexing from is:ie */

  double lx1 = 1.0; //these values are inclusive [x1_i, x1_f]
  double lx2 = 1.0; 
  double x1_i = 0.0;
  double x2_i = 0.0;

  double dx2 = lx2/(nx2_r-1);
  double x2_f = x2_i + lx2;

  double dx1 = lx1/(nx1_r-1);   
  double x1_f = x1_i + lx1;
  
  /*Cell centered (zonal) values of computational coordinate position */
  double *x1 = (double *) malloc(sizeof(double)*nx1); 
  double *x2 = (double *) malloc(sizeof(double)*nx2);
  x1[is] = x1_i;
  x2[js] = x2_i;
  for(i=is+1; i<ie; i++){ 
    x1[i] = x1[i-1] + dx1;
  }
  for(i=js+1; i<je; i++){
    x2[i] = x2[i-1] + dx2;
  } 
    
  /* Mesh edge (nodal) values of computational coordinate position */
  /* Nodal quantities refer to the left state of the cell, i.e. i-th index is i-1/2 boundary */
  double *x1_b = (double *) malloc(sizeof(double)*(nx1+1)); 
  double *x2_b = (double *) malloc(sizeof(double)*(nx2+1));
  x1_b[is] = x1_i - dx1/2;
  x2_b[js] = x2_i - dx2/2;
  for(i=is+1; i<=ie; i++){ 
    x1_b[i] = x1_b[i-1] + dx1;
  }
  for(i=js+1; i<=je; i++){
    x2_b[i] = x2_b[i-1] + dx2;
  } 

  /*Cell centered (zonal) values of physical coordinate position */
  //These must be 2D arrays since coordinate transformation is not diagonal
  //indexed by x1 columns and x2 rows.

  /*CONVENTION: flip the dimension ordering in all multiD arrays due to row major ordering of C */
  double *dataX = (double *) malloc(sizeof(double)*nx1*nx2);
  double *dataY = (double *) malloc(sizeof(double)*nx1*nx2);
  double **x = allocate_2D_contiguous(dataX,nx2,nx1); 
  double **y = allocate_2D_contiguous(dataY,nx2,nx1); 

  /*precompute coordinate mappings */
  for(j=js; j<je; j++){
    for(i=is; i<ie; i++){
      x[j][i] = x1[i];
      y[j][i] = x2[j]; 
    }
  }

  /*Mesh edge (nodal) values of physical coordinate position */
  //This is a relic of non-Carteisan code
  double *dataXb = (double *) malloc(sizeof(double)*(nx1+1)*(nx2+1));
  double *dataYb = (double *) malloc(sizeof(double)*(nx1+1)*(nx2+1));
  double **x_b = allocate_2D_contiguous(dataXb,nx2+1,nx1+1); 
  double **y_b = allocate_2D_contiguous(dataYb,nx2+1,nx1+1); 
  for(j=js; j<=je; j++){
    for(i=is; i<=ie; i++){
      x_b[j][i] = x1_b[i];
      y_b[j][i] = x2_b[j]; 
    }
  }
 
  /*----------------------------------------*/
  /*         Initialize Angular Mesh        */
  /*----------------------------------------*/

  /*Discretize photon directions */
  double *xa1 = (double *) malloc(sizeof(double)*nxa1); 
  double *xa1_b = (double *) malloc(sizeof(double)*(nxa1+1)); 
  double *datamu = (double *) malloc(sizeof(double)*nxa1*3); 
  double *datamu_b = (double *) malloc(sizeof(double )*3*(nxa1+1));//identify with left and bottom boundaries 
  double *pw = (double *) malloc(sizeof(double)*nxa1); 
  double **mu = allocate_2D_contiguous(datamu,nxa1,3); //CONVENTION: xa1 is the first dimension, component is second dim 
  double **mu_b = allocate_2D_contiguous(datamu_b,nxa1+1,3); 

  if (xa1_uniform){
    nxa1 = uniform_angles2D(nxa1, phi_z, pw, mu, mu_b, xa1,xa1_b); 
  }
  else{
    bruls_angles2D(N_bruls, phi_z, pw, mu, mu_b, xa1,xa1_b); 
  }

  for (i=0; i<nxa1-1; i++)
    dxa1[i] = xa1_b[i+1] - xa1_b[i];
  //last boundary should be 2pi 
  dxa1[nxa1-1] = 2*M_PI - xa1_b[nxa1-1];

  /*----------------------------------------*/
  /*    Precompute advection parameters     */
  /*----------------------------------------*/
  printf("Precomputing advection parameters:\n");
  /* Compute phase volumes of all cells */
  double *datakappa= (double *) malloc(sizeof(double)*nx1*nx2*nxa1);
  double ***vol = allocate_3D_contiguous(datakappa,nxa1+1,nx2,nx1); 
  for(k=ks; k<=ke; k++){
    for(j=js-1; j<=je; j++){
      for(i=is-1; i<=ie; i++){	
	vol[k][j][i] = dx1*dx2; 
	if (k==ke)
	  vol[k][j][i] *= dxa1[ks];
	else
	  vol[k][j][i] *= dxa1[k];	
      }
    }
  } 

  /* Compute average normal edge velocities (or approximations) */
  // Again, these are nodal values, so the index i refers to i-1/2 value
  double *dataU = malloc(sizeof(double)*nxa1*nx2*nx1);
  double *dataV = malloc(sizeof(double)*nxa1*nx2*nx1);

  double ux,vy,temp; 
  double ***U = allocate_3D_contiguous(dataU,nxa1+1,nx2,nx1);
  double ***V = allocate_3D_contiguous(dataV,nxa1+1,nx2,nx1);
  double ***source = allocate_3D_contiguous(dataV,nxa1+1,nx2,nx1);

  for(k=ks; k<=ke; k++){  
    if (k==ke) //need this to compute exact average edge velocity integral
      next = ks+1;
    else
      next = k+1;
    for(j=js-1; j<=je; j++){       
      for(i=is-1; i<=ie; i++){
	/* Midpoint approximations to edge velocities  */
	if (k==ke){
	  U[k][j][i] = mu[ks][0]; // cos(xa1_b[next]); 
	  V[k][j][i] = mu[ks][1]; //mu[next][0]*sin(dx2/2) + mu[next][1]*cos(dx2/2); //sin(xa1[next] +dx2/2); //
	}
	else{
	  U[k][j][i] = mu[k][0]; // cos(xa1_b[next]); 
	  V[k][j][i] = mu[k][1]; //mu[next][0]*sin(dx2/2) + mu[next][1]*cos(dx2/2); //sin(xa1[next] +dx2/2); //
	}	  
	/* Compute exact average edge velocities */
	// ke is a ghost cell necessary for proper Visit output
	//dxa1[ke] not defined
	//xa1_b[ks] = xa1_b[ke], and we want the velocity to match the first real cell 
#ifdef V_EXACT
	U[k][j][i] = (sin(xa1_b[next]) - sin(xa1_b[k]))/dxa1[ks];
	V[k][j][i] = (-cos(xa1_b[next]) + cos(xa1_b[k]))/dxa1[ks];	   
#endif
	source[k][j][i] = 0.0;
      }
    }
  }
  //  printf("xa1_b[ks] = %lf xa1_b[ke] = %lf dxa1[ke] = %lf \n",xa1_b[ks],xa1_b[ke],dxa1[ke]); 
  /* Print out all velocities */
  /*  for(k=ks; k<=ke; k++){  
    for(j=js-1; j<=je; j++){       
      for(i=is-1; i<=ie; i++){
	printf("(U,V)[%d,%d,%d] = (%lf, %lf) \n",k,j,i,U[k][j][i],V[k][j][i]); 
	}}}  */

  /* Use numerical integration to compute solid angle projections for #ifdef SOLID_ANGLE case */
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000); //max number of subintervals per integral
  double result, error;
  double alpha;  //no params called by function

  gsl_function F;
  F.params = &alpha;

  //average edge velocity of projection of solid angle cone with positive projection on cell
  // boundary for uniform carteisn discretization
  double *solid_U_plus = malloc(sizeof(double)*nxa1); 
  double *solid_U_minus = malloc(sizeof(double)*nxa1); 
  double *solid_V_plus = malloc(sizeof(double)*nxa1); 
  double *solid_V_minus = malloc(sizeof(double)*nxa1); 

  for (k=0; k<nxa1; k++){
    F.function = &max_cosine;
    gsl_integration_qags (&F, xa1_b[k], xa1_b[k]+dxa1[k], 0, 1e-7, 1000,w, &result, &error); 
    solid_U_plus[k] = result/dxa1[k]; 
    F.function = &min_cosine;
    gsl_integration_qags (&F, xa1_b[k], xa1_b[k]+dxa1[k], 0, 1e-7, 1000,w, &result, &error); 
    solid_U_minus[k] = result/dxa1[k]; 
    F.function = &max_sine;
    gsl_integration_qags (&F, xa1_b[k], xa1_b[k]+dxa1[k], 0, 1e-7, 1000,w, &result, &error); 
    solid_V_plus[k] = result/dxa1[k]; 
    F.function = &min_sine;
    gsl_integration_qags (&F, xa1_b[k], xa1_b[k]+dxa1[k], 0, 1e-7, 1000,w, &result, &error); 
    solid_V_minus[k] = result/dxa1[k];
  }

  /* Check CFL condition, reset timestep */
  double *datacfl = malloc(sizeof(double)*nxa1*nx2*nx1); 
  double ***cfl_array =  allocate_3D_contiguous(datacfl,nxa1,nx2,nx1);

  for(k=0; k<nxa1; k++){ //based on edge velocities 
    for(j=0; j<nx2; j++){
      for(i=0; i<nx1; i++){
	if (i >=is && i< ie && j >=js && j <je)
	  cfl_array[k][j][i] = fabs(U[k][j][i])*dt/dx1 + fabs(V[k][j][i])*dt/(dx2); 
	else
	  cfl_array[k][j][i] =0.0; 
	datacfl[index] = cfl_array[k][j][i];
	index++;
      }
    }
  }

  //find maximum CFL value in domain
  double max_cfl = find_max_double(datacfl,nxa1*nx1*nx2); 
  index = 0;
  printf("Largest CFL number = %lf\n",max_cfl); 
  if (max_cfl > CFL || AUTO_TIMESTEP){ //reset timestep if needed
    if (max_cfl ==0){
      printf("Max CFL is 0.0!\n"); //dont divide by 0
      exit(1);
    }
    printf("Original dt = %lf\n",dt);
    dt = CFL*dt/max_cfl; 
    printf("Auto-timestep dt = %lf\n",dt);
    for(k=0; k<nxa1; k++){ //ERROR IN ADVECTION.C: ARRAY STARTS AT 1
      for(j=0; j<nx2; j++){ //POSSIBLE ERROR IN POLAR COORDINATES SOLVER-- CHECK BRACES WITH FOR LOOPS AND IF STATEMENTS
	for(i=0; i<nx1; i++){ 
	  if (i >=is && i< ie && j >=js && j <je)
	    cfl_array[k][j][i] = fabs(U[k][j][i])*dt/dx1 + fabs(V[k][j][i])*dt/(dx2);
	  else
	    cfl_array[k][j][i] =0.0; 	
	  datacfl[index] = cfl_array[k][j][i];
	  index++;
	}
      } 
    }
  }
  /* Confirm correct CFL number */
  max_cfl = find_max_double(datacfl,nxa1*nx1*nx2);
  printf("Largest CFL number = %.5e\n",max_cfl); 
  
  /*Conserved variable on the computational coordinate mesh*/
  /* Solid angle averaged, spatial averaged intensity */
  double *dataI =  malloc(sizeof(double)*nx1*nx2*nxa1);
  double ***I = allocate_3D_contiguous(dataI,nxa1,nx2,nx1); 

  /* Net fluxes in each dimension at each timestep */
  double U_plus,U_minus,V_plus,V_minus; // X_plus = max(0.0, X), e.g. //ERROR: these are deprecated
  double *dataFlux = malloc(sizeof(double)*nx1*nx2*nxa1);
  double ***net_flux = allocate_3D_contiguous(dataFlux,nxa1,nx2,nx1); 

  /*----------------------------------------*/
  /*         Configure output               */
  /*----------------------------------------*/
  printf("Setting up VTK output....\n");
  /* Using Visit VTK writer */
  char filename[20];
  /*CONVENTION: For the 2D1V solver, the angular dimension is treated as third spatial dimension.
    May want to create separate variables in the future */
  int dims[] = {nx1_r+1, nx2_r+1, nxa1+1}; //dont output ghost cells. nodal variables have extra edge point
  int nvars = 3;
  int vardims[] = {1, 1, 3}; //I is a scalar, J is a scalar, edge velocity is a 3-vector 
  int centering[] = {0, 0, 1}; // I,J are cell centered (zonal), velocity is defined at edges (nodal)
  const char *varnames[] = {"I","J", "vel"};
  /* Curvilinear mesh points stored x0,y0,z0,x1,y1,z1,...*/
  //An array of size nI*nJ*nK*3 . These points are nodal, not zonal
  float *pts = (float *) malloc(sizeof(float)*(nx1_r+1)*(nx2_r+1)*(nxa1+1)*3); //check angular dimension size
  //The array should be layed out as (pt(i=0,j=0,k=0), pt(i=1,j=0,k=0), ...
  //pt(i=nI-1,j=0,k=0), pt(i=0,j=1,k=0), ...).
  index=0; 
  for(k=ks; k<=ke; k++){
    for(j=js; j<=je; j++){
      for(i=is; i<=ie; i++){
	pts[index] = x_b[j][i];
	pts[++index] = y_b[j][i];
	pts[++index] = xa1_b[k]; 
	index++;
      }
    }
  }
  
  /* pack U,V,W into a vector */
  float *edge_vel = (float *) malloc(sizeof(float)*(nx1_r+1)*(nx2_r+1)*(nxa1+1)*3); //An array of size nI*nJ*nK*3 
  index=0; 
  for(k=ks; k<=ke; k++){
    for(j=js; j<=je; j++){ //ERROR: j=je ghost cells are messed up since U,V,W arent initialized that far
      for(i=is; i<=ie; i++){
	edge_vel[index] = U[k][j][i];
	edge_vel[++index] = V[k][j][i]; 
	edge_vel[++index] = 0.0; 
	index++;
      }
    }
  } 

  //  vars       An array of variables.  The size of vars should be nvars.
  //                 The size of vars[i] should be npts*vardim[i].
  float *realI, *realJ; //these flattened 1D vectors exclude spatial ghost cells
  realI =(float *) malloc(sizeof(float)*nx1_r*nx2_r*nxa1);//ERROR IN ADVECTION. SIZEOF(DOUBLE)
  realJ =(float *) malloc(sizeof(float)*nx1_r*nx2_r*nxa1); //unfortunately, if it lives on the same 3D mesh, must be duplicated at each z/nxa1 height
  float *vars[] = {(float *) realI,(float *) realJ, (float *)edge_vel};

  /*----------------------------------------*/
  /*  Compute and output analytic solution  */
  /*----------------------------------------*/
#ifdef ANALYTIC_SOLUTION
  printf("Writing analytic solution....\n"); 
  float *realI_analytic, *realJ_analytic; 
  realI_analytic =(float *) malloc(sizeof(float)*nx1_r*nx2_r*nxa1);
  realJ_analytic =(float *) malloc(sizeof(float)*nx1_r*nx2_r*nxa1);
 
  index=0;
  double x_source, y_source, slope_c,slope_min, slope_max; 
  double dist; //distance from source with radius dx/2

  /* First, plot the analytic solution for the diagonal search light (with one source point and \Delta \xi = \pi / 2 */
  //relevant indices of problem
  int source_i = is + nx1_r/2;
  int source_j = js + nx2_r/2;
  double source_x = x1[source_i];
  double source_y = x2[source_j];
  double source_radius = dx1/2; 

  /* Simplified analytic solution: zero out all cells (since search algorithm increments intensity bin) */
  for (k=0; k<ke; k++){
    for (j=js; j<je; j++){
      for (i=is; i<ie; i++){
	I[k][j][i] = 0.0; 
      }
    }
  }

  //source in bottom left corner
  /*  source_x = x1_b[source_i];
      source_y = x2_b[source_j]; */

  //centered source-- this breaks the diagonal beam case, since some rays correspond to angles 3pi/2<xi<2pi. Use in isotropic case
  source_x = x1[source_i];
  source_y = x2[source_j];

  int subsamples = 10; //= number of points inside a cell in 1D
  double dx_s = dx1/(subsamples); 
  double dy_s = dx2/(subsamples); 
  double rec_angle; 
  double rec_x, rec_y; 
  int debug_index = 0; 
  index =0; 
  for (j=js; j<je; j++){
    for (i=is; i<ie; i++){ 
      if (i==source_i && j==source_j)
	for (k=source_k_start; k<=source_k_end; k++){
	  I[k][j][i] = 1.0;
	}      
      else{
	//	if (i >= source_i && j >= source_j){ //ONLY USE IN DIAGONAL BEAM TEST: restricts to first quadrant. Comment out lower brace
	  //Simplified approximation to analytic solution: subsample each cell spatially 
	  for (l=0; l<subsamples; l++){
	    for (n=0; n<subsamples; n++){ //is this the correct way to calculate subsample cell centers?
	      rec_x = (x1_b[i] + dx_s/2) // first subsample x position
		+ l*dx_s; 
	      rec_y = (x2_b[j] + dy_s/2) // first subsample y position
		+ n*dy_s; 
	      dist = sqrt((rec_x - source_x)*(rec_x - source_x) + (rec_y - source_y)*(rec_y - source_y));
	      rec_angle = atan2(rec_y-source_y,rec_x-source_x);       //4 quadrant arctangent returns values -pi, pi
	      rec_angle = (rec_angle > 0 ? rec_angle : (2*M_PI + rec_angle)) ;       // if angle is less than 0 radians, add 2pi radians
	      //		printf("distance = %lf angle = %lf\n",dist,rec_angle);
	      debug_index = 0; 
	      //place intensity in corresponding solid angle 
	      for (k=source_k_start; k<=source_k_end; k++){
		//		  printf("xa1_b[%d] = %lf, xa1_b[%d] = %lf\n",k,xa1_b[k],k+1,xa1_b[k+1]); 
		if (rec_angle >= xa1_b[k] && rec_angle <= xa1_b[k+1]){ // if the angle of the sample point is equal to a boundary, split the intensity in half.
		  debug_index = 1; 
		  I[k][j][i] += (dx1/2)*(1.0/dist);
		  // debug first row in first quadrant-- angle should be near zero radians
		  /*		    if (i==source_i+1 && j==source_j){
				    printf("cell (%d,%d) subsample (%d, %d) found angular bin k=%d \n",i,j,l,n,k); 
				    printf("x,y = (%lf, %lf) distance = %lf angle = %lf\n",rec_x,rec_y,dist,rec_angle);
				    }*/
		}
	      }
	      if (debug_index == 0){ //ensure that every subsample is being placed in some solid angular bin
		printf("angular bin not found!\n"); 
		return(1); 
	      }
	    }
	  }
	  //im not sure that our phase volume average falls off as 1/r for a uniform discretized cell... actually I think it does. 
	  for (k=0; k<ke; k++){//renormalize the phase space cell
	    I[k][j][i] /= (subsamples*subsamples)*dxa1[k]; //need to weight by dxa1 otherwise falloff is too fast
	  }
	  //} //ONLY USE IN DIAGONAL TEST CASE: comment out for isotropic test case
      }
    }
  }
  for (k=ks; k<ke; k++){
    for (j=js; j<je; j++){
      for (i=is; i<ie; i++){
	realI_analytic[index] = (float) I[k][j][i];
	index++;
      }
    }
  }

  /* Compute zeroth moment of analytic radiation field */
  indexJ=0;
  for (j=0; j<nx2_r; j++){   //zero out J array:
    for (i=0; i<nx1_r; i++){
      realJ_analytic[indexJ] = 0;
      indexJ++;
    }
  }
  index = 0; 
  for (k=ks; k<ke; k++){
    indexJ=0; 
    for (j=js; j<je; j++){
      for (i=is; i<ie; i++){
	realJ_analytic[indexJ] += (float) realI_analytic[index]*pw[k];
	indexJ++;
	index++;
      }
    }
  }
  
//Debug change in analytic J with N_\xi
/*    indexJ=0; 
    for (j=js; j<je; j++){
      for (i=is; i<ie; i++){
	if (i==source_i+1 && j==source_j){
	  printf("J[%d][%d] = %lf\n",j,i,realJ[indexJ]);
	  for (k=0; k<ke; k++){
	    printf("I[%d][%d][%d] = %lf pw = %lf xi: [%lf, %lf] \n",k,j,i,I[k][j][i],pw[k],xa1_b[k],xa1_b[k+1]);
	  }}	  
	indexJ++;
      }
    }
*/
  
  /* Copy the analytic radiation energy field from first angular bin to all angular bins */
  indexJ=nx1_r*nx2_r; 
  for (k=ks+1; k<ke; k++){
    for (j=0; j<nx2_r; j++){
      for (i=0; i<nx1_r; i++){
	realJ_analytic[indexJ] = realJ_analytic[j*nx1_r+i]; 
	indexJ++;
      }
    }
  }

  sprintf(filename,"analytic-rte.vtk"); 

  //kyle: new analytic solution: temporarily change pointers before doing actual computation:
  vars[0] = realI_analytic;
  vars[1] = realJ_analytic;
  write_curvilinear_mesh(filename,1,dims, pts, nvars,vardims, centering, varnames, vars);
  //  return(0); //uncomment if you just want analytic solution
#endif //ANALYTIC_SOLUTION

  /* Reset output pointers to numerical solution */
  vars[0] = realI; 
  vars[1] = realJ;

  /* Compute initial condition */
  for(k=ks; k<ke; k++){
    for(j=js; j<je; j++){
      for(i=is; i<ie; i++){
	I[k][j][i] = initial_condition(x[j][i],y[j][i],xa1[k]); 
      }
    }
  }

  /* Output initial condition */
  index=0; 
  indexJ=0;
  for (j=0; j<nx2_r; j++){   //zero out J array:
    for (i=0; i<nx1_r; i++){
      realJ[indexJ] = 0;
      indexJ++;
    }
  }
  for (k=ks; k<ke; k++){
    indexJ=0; 
    for (j=js; j<je; j++){
      for (i=is; i<ie; i++){
	realI[index] = (float) I[k][j][i];
	realJ[indexJ] += (float) I[k][j][i]*pw[k];
	index++;
	indexJ++;
      }
    }
  }

  /* Copy the radiation energy field from first angular bin to all angular bins */
  indexJ=nx1_r*nx2_r; 
  for (k=ks+1; k<ke; k++){
    for (j=0; j<nx2_r; j++){
      for (i=0; i<nx1_r; i++){
	realJ[indexJ] = realJ[j*nx1_r+i]; 
	indexJ++;
      }
    }
  }
  sprintf(filename,"rte-000.vtk"); 
  write_curvilinear_mesh(filename,1,dims, pts, nvars,vardims, centering, varnames, vars);

#ifdef FACE_INTERPOLATE
  //should reuse this structures, only keeping final flux
  double ***da1 = allocate_3D_contiguous(dataFlux,nxa1,nx2,nx1);
  double ***da2 = allocate_3D_contiguous(dataFlux,nxa1,nx2,nx1);
  double ***a1_l = allocate_3D_contiguous(dataFlux,nxa1,nx2,nx1); //temporary left x-edge face extrapolated values
  double ***a1_r = allocate_3D_contiguous(dataFlux,nxa1,nx2,nx1); //temporary left x-edge face extrapolated values
  double ***a2_l = allocate_3D_contiguous(dataFlux,nxa1,nx2,nx1); //temporary left x-edge face extrapolated values
  double ***a2_r = allocate_3D_contiguous(dataFlux,nxa1,nx2,nx1); //temporary left x-edge face extrapolated values
  double a6; 
  double ***xFace_I = allocate_3D_contiguous(dataFlux,nxa1,nx2,nx1); 
  double ***yFace_I = allocate_3D_contiguous(dataFlux,nxa1,nx2,nx1); 
  double ***xFace_flux = allocate_3D_contiguous(dataFlux,nxa1,nx2,nx1); 
  double ***yFace_flux = allocate_3D_contiguous(dataFlux,nxa1,nx2,nx1); 
#endif

  /*-----------------------*/
  /* Main timestepping loop */
  /*-----------------------*/
  printf("dt=%lf dx1=%lf dx2=%lf dxa1=%lf \n",dt,dx1,dx2,dxa1[0]); 
//  printf("dt=%.5e dx1=%.5e dx2=%.5e dxa1=%.5e \n",dt,dx1,dx2,dxa1[0]); 
  for (n=1; n<nsteps; n++){
    /*Spatial boundary conditions */
    //bcs are specified along a computational coord direction, but are a function of the physical coordinate of adjacent "real cells"
    for(l=0;l<num_ghost; l++){
      for(k=ks;k<ke; k++){
	for (j=js; j<je; j++){
	  if(X1_PERIODIC){
	    I[k][j][l] = I[k][j][ie-1-l];
	    I[k][j][nx1-1-l] = I[k][j][is+l];
	  }
	  else{
	    I[k][j][l] = bc_x1i(x[j][is],y[j][is],xa1[k],n*dt);
	    I[k][j][nx1-1-l] = bc_x1f(x[j][is],y[j][is],xa1[k],n*dt);
	  }
	}	      
	for (i=is; i<ie; i++){
	  if(X2_PERIODIC){
	    I[k][l][i] = I[k][je-1-l][i];
	    I[k][nx2-1-l][i] = I[k][js+l][i];
	  }
	  else{ 
	    I[k][l][i] = bc_x2i(x[js][i],y[js][i],xa1[k],n*dt);
	    I[k][nx2-1-l][i] = bc_x2f(x[je-1][i],y[je-1][i],xa1[k],n*dt);
	  }  
	} 
      }
    }
    
    /* set any fixed points inside domain */
    double fixed_I; 
    for (k=ks; k<ke; k++){
      for (j=js; j<je; j++){
	for (i=is; i<ie; i++){
#ifdef TEST5
	  if (i == is + nx1_r/2 && j == js + nx2_r/2 && k >= source_k_start && k <= source_k_end){ //integer division
	    I[k][j][i] = 1.0;	   	
	  }
#endif 
	  fixed_I = bc_interior(x1[i], x2[j], xa1[k], n*dt);
	  if (fixed_I != 0.0)
	    I[k][j][i] = fixed_I; 
	}
      }
    }
    
    /* SOLID_ANGLE method requires fluxes to manually be cleared */
    for (k=ks; k<ke; k++){
      for (j=js; j<je; j++){
	for (i=is; i<ie; i++){
	  net_flux[k][j][i] =0.0; 
	}
      }
    }

#ifdef FACE_INTERPOLATE
    /* At each timestep, compute the face averages of the intensity with fourth order central, linear interpolation limited by PPM */
    //compute average linear slopes, need to do for first ghost cells
    for (k=ks; k<ke; k++){
      for (j=js-1; j<=je; j++){
	for (i=is-1; i<=ie; i++){	 
	  da1[k][j][i] = 0.5*(I[k][j][i+1] - I[k][j][i]);
	  da2[k][j][i] = 0.5*(I[k][j+1][i] - I[k][j][i]);
	}
      }
    }

    // monotonize slopes
    for (k=ks; k<ke; k++){
      for (j=js-1; j<=je; j++){
	for (i=is-1; i<=ie; i++){	 
	  if ((I[k][j][i+1] - I[k][j][i])*(I[k][j][i] - I[k][j][i-1]) > 0){
	    temp = fmin(fabs(da1[k][j][i]),2*fabs(I[k][j][i] - I[k][j][i-1]));
	    da1[k][j][i] = fmin(temp,fabs(I[k][j][i+1] - I[k][j][i])); 
	  }
	  else 
	    da1[k][j][i] = 0.0;
	  if ((I[k][j+1][i] - I[k][j][i])*(I[k][j][i] - I[k][j-1][i]) > 0){
	    temp = fmin(fabs(da2[k][j][i]),2*fabs(I[k][j][i] - I[k][j-1][i]));
	    da2[k][j][i] = fmin(temp,fabs(I[k][j+1][i] - I[k][j][i])); 
	  }
	  else 
	    da2[k][j][i] = 0.0;
	}
      }
    }

    //compute l/r states
    for (k=ks; k<ke; k++){
      for (j=js; j<je; j++){
	for (i=is; i<ie; i++){	  //a1_l(i) may not equal a1_r(i-1) after limiting
	  //these expressions reference da1 in ghost cells, which are not initialized?
	  a1_l[k][j][i] = I[k][j][i-1] +0.5*(I[k][j][i] - I[k][j][i-1]) + (da1[k][j][i-1] - da1[k][j][i])/6; 
	  a1_r[k][j][i] = I[k][j][i] +0.5*(I[k][j][i+1] - I[k][j][i]) + (da1[k][j][i] - da1[k][j][i+1])/6; 

	  a2_l[k][j][i] = I[k][j-1][i] +0.5*(I[k][j][i] - I[k][j-1][i]) + (da2[k][j-1][i] - da2[k][j][i])/6; 
	  a2_r[k][j][i] = I[k][j][i] +0.5*(I[k][j+1][i] - I[k][j][i]) + (da2[k][j][i] - da2[k][j+1][i])/6; 
	  //when solution is smooth on a uniform grid, these are equal to:
	  /*
	  a1_l[k][j][i] = 7.0/12.0 * (I[k][j][i-1] + I[k][j][i]) - 1.0/12.0*(I[k][j][i-2] +I[k][j][i+1]);
	  a1_r[k][j][i] = 7.0/12.0 * (I[k][j][i] + I[k][j][i+1]) - 1.0/12.0*(I[k][j][i-1] +I[k][j][i+2]);
	  a2_l[k][j][i] = 7.0/12.0 * (I[k][j-1][i] + I[k][j][i]) - 1.0/12.0*(I[k][j-2][i] +I[k][j+1][i]);
	  a2_l[k][j][i] = 7.0/12.0 * (I[k][j][i] + I[k][j+1][i]) - 1.0/12.0*(I[k][j-1][i] +I[k][j+2][i]); */
	}
      }
    }						

    // monotonize l/r states
    for (k=ks; k<ke; k++){
      for (j=js; j<je; j++){
	for (i=is; i<ie; i++){       
	  //enforce monotonicity in x. sekora collela eq 14, ppm paper eq 1.10
	  if ((a1_r[k][j][i] - I[k][j][i])*(I[k][j][i] - a1_l[k][j][i]) <= 0){ //sekora and original ppm article seem to disagree here
	    a1_l[k][j][i] = I[k][j][i]; 
	    a1_r[k][j][i] = I[k][j][i]; 
	  }
	  else if ((a1_r[k][j][i] - a1_l[k][j][i])*(I[k][j][i] - 0.5*(a1_r[k][j][i] + a1_l[k][j][i])) > (a1_r[k][j][i] - a1_l[k][j][i])*(a1_r[k][j][i] - a1_l[k][j][i])/6){
	    a1_l[k][j][i] = 3*I[k][j][i] - 2*a1_r[k][j][i]; 
	    //	    if (k==0 && ((i==18 && j==19) || (i==19 && j==18)))
	           // printf("overshoot on left a1\n");
	    // this is sekora eq 14
	    //(fabs(a1_l[k][j][i] - I[k][j][i]) >= 2*fabs(a1_r[k][j][i] - I[k][j][i])){ //overshoot of interpolant by a_L
	    //	  a1_l[k][j][i] = I[k][j][i] - 2*(a1_r[k][j][i] - I[k][j][i]); //agrees with ppm
	  }
	  else if ((a1_r[k][j][i] - a1_l[k][j][i])*(I[k][j][i] - 0.5*(a1_r[k][j][i] + a1_l[k][j][i])) < -(a1_r[k][j][i] - a1_l[k][j][i])*(a1_r[k][j][i] - a1_l[k][j][i])/6){
	    a1_r[k][j][i] = 3*I[k][j][i] - 2*a1_l[k][j][i]; 	   
	  }

	  //enforce monotonicity in y
	  if ((a2_r[k][j][i] - I[k][j][i])*(I[k][j][i] - a2_l[k][j][i]) <= 0){ //sekora and original ppm article seem to disagree here
	    a2_l[k][j][i] = I[k][j][i]; 
	    a2_r[k][j][i] = I[k][j][i]; 
	  }
	  else if ((a2_r[k][j][i] - a2_l[k][j][i])*(I[k][j][i] - 0.5*(a2_r[k][j][i] + a2_l[k][j][i])) > (a2_r[k][j][i] - a2_l[k][j][i])*(a2_r[k][j][i] - a2_l[k][j][i])/6){
	    a2_l[k][j][i] = 3*I[k][j][i] - 2*a2_r[k][j][i]; 
	  }
	  else if ((a2_r[k][j][i] - a2_l[k][j][i])*(I[k][j][i] - 0.5*(a2_r[k][j][i] + a2_l[k][j][i])) < -(a2_r[k][j][i] - a2_l[k][j][i])*(a2_r[k][j][i] - a2_l[k][j][i])/6){
	    a2_r[k][j][i] = 3*I[k][j][i] - 2*a2_l[k][j][i]; 
	  }
	}
      }
    }

    //BIG QUESTION: in the banks + hittinger paper, do we compute separate edge fluxes for L/R states?? or should they be identical. How do we reconcile that with the possibility that I_L/R may not be equal on a cell edge? ANSWER: Solving the riemann problem selects the L/R state of the primitive variable. No L/R fluxes.... only one flux per boundary in order to satisfy conservation

    //so, let's try this. Following McCorquodale and Colella, from L/R states, solve riemann problem and pick one edge variable
    for (k=ks; k<ke; k++){
      for (j=js; j<je; j++){
	for (i=is; i<ie; i++){
	  if (U[k][j][i] >= 0.0)
	    xFace_I[k][j][i] = a1_r[k][j][i-1]; 
	  else 
	    xFace_I[k][j][i] = a1_l[k][j][i]; 
	  if (V[k][j][i] >= 0.0)
	    yFace_I[k][j][i] = a2_r[k][j-1][i]; 
	  else 
	    yFace_I[k][j][i] = a2_l[k][j][i]; 
	}
      }
    }

    for (k=ks; k<ke; k++){
      for (j=js; j<je; j++){
	for (i=is; i<ie; i++){
	  if (k==ke-1)
	    next = ks; 
	  else 
	    next = k+1; 
	  if (k==ks)
	    prev = ke-1;
	  else
	    prev = k-1; 
	  //temporarily turn off the transverse terms 
	  //old L/R flux formulation
	  //	  xFace_flux_l[k][j][i] = U[k][j][i]*a1_l[k][j][i];// + 1.0/48.0 *(a1_l[next][j][i] - a2_l[prev][j][i])*(U[next][j][i] - U[prev][j][i]); 
	  //xFace_flux_r[k][j][i] = U[k][j][i]*a1_r[k][j][i];// + 1.0/48.0 *(a1_r[next][j][i] - a2_r[prev][j][i])*(U[next][j][i] - U[prev][j][i]); 
	  /*
	  if (xFace_flux_l[k][j][i] != 0.0)
	    printf("%d,%d,%d x_flux_l = %lf (I_x)_L = %lf \n",k,j,i,xFace_flux_l[k][j][i],a1_l[k][j][i]);
	  if (yFace_flux_l[k][j][i] != 0.0)
	    printf("%d,%d,%d y_flux_l = %lf (I_y)_L = %lf \n",k,j,i,yFace_flux_l[k][j][i],a2_l[k][j][i]);
	  if (xFace_flux_r[k][j][i] != 0.0)
	    printf("%d,%d,%d x_flux_r = %lf (I_x)_R = %lf \n",k,j,i,xFace_flux_r[k][j][i],a1_r[k][j][i]);
	  if (yFace_flux_r[k][j][i] != 0.0)
	  printf("%d,%d,%d y_flux_r = %lf (I_y)_R = %lf \n",k,j,i,yFace_flux_r[k][j][i],a2_r[k][j][i]);   
	  */

	  //original PPM formulation for 1D advection:
	  /*	  if (U[k][j][i] >= 0.0){
	    //center at i-1
	    a6= 6*(I[k][j][i-1] - 0.5*(a1_l[k][j][i-1] + a1_r[k][j][i-1])); 
	    xFace_flux[k][j][i] = U[k][j][i]*dt/dx1*(xFace_I[k][j][i] - U[k][j][i]*dt/dx1 *((a1_r[k][j][i-1] - a1_l[k][j][i-1]) - (1-2.0/3.0 * U[k][j][i]*dt/dx1)*a6));
	  }
	  else{
	    //center at i
	    a6= 6*(I[k][j][i] - 0.5*(a1_l[k][j][i] + a1_r[k][j][i])); 
	    xFace_flux[k][j][i] = U[k][j][i]*dt/dx1*(xFace_I[k][j][i] - U[k][j][i]*dt/dx1 *((a1_r[k][j][i] - a1_l[k][j][i]) + (1+2.0/3.0 * U[k][j][i]*dt/dx1)*a6));
	  }
	  if (V[k][j][i] >= 0.0){
	    //center at j-1
	    a6= 6*(I[k][j-1][i] - 0.5*(a2_l[k][j-1][i] + a2_r[k][j-1][i])); 
	    yFace_flux[k][j][i] = V[k][j][i]*dt/dx2*(yFace_I[k][j][i] - V[k][j][i]*dt/dx2 *((a2_r[k][j-1][i] - a2_l[k][j-1][i]) - (1-2.0/3.0 * V[k][j][i]*dt/dx2)*a6));
	  }
	  else{
	    //center at j
	    a6= 6*(I[k][j][i] - 0.5*(a2_l[k][j][i] + a2_r[k][j][i])); 
	    yFace_flux[k][j][i] = V[k][j][i]*dt/dx2*(yFace_I[k][j][i] - V[k][j][i]*dt/dx2 *((a2_r[k][j][i] - a2_l[k][j][i]) + (1+2.0/3.0 * V[k][j][i]*dt/dx2)*a6));
	    } */

	  //new consistent edge flux formulation
	  xFace_flux[k][j][i] = U[k][j][i]*xFace_I[k][j][i] + 1.0/48.0 *(xFace_I[next][j][i] - xFace_I[prev][j][i])*(U[next][j][i] - U[prev][j][i]); 
	  yFace_flux[k][j][i] = V[k][j][i]*yFace_I[k][j][i] + 1.0/48.0 *(yFace_I[next][j][i] - yFace_I[prev][j][i])*(V[next][j][i] - V[prev][j][i]);  

	  //debug statements
	  /*	  if (xFace_flux[k][j][i] != 0.0)
	    printf("%d,%d,%d x_flux = %lf x_I = %lf (I_x-1)_R = %lf (I_x)_L = %lf \n",k,j,i,xFace_flux[k][j][i],xFace_I[k][j][i], a1_r[k][j][i-1],a1_l[k][j][i]);
	  if (yFace_flux[k][j][i] != 0.0)
	  printf("%d,%d,%d y_flux = %lf y_I = %lf (I_y-1)_R = %lf (I_y)_L = %lf \n",k,j,i,yFace_flux[k][j][i],yFace_I[k][j][i],a2_r[k][j-1][i],a2_l[k][j][i]);	 */
	}
      }
    }

#endif

    double flux_limiter =0.0; 
    double flux_limiter_x1_l,flux_limiter_x1_r,flux_limiter_x2_l,flux_limiter_x2_r,flux_limiter_xa1_l,flux_limiter_xa1_r; 
    double *imu = (double *) malloc(sizeof(double)*3); //manually copy array for computing slope limiters
    double *pr = (double *) malloc(sizeof(double)*3); //manually copy array for computing slope limiters
    /* Donor cell upwinding */
    for (k=ks; k<ke; k++){
      for (j=js; j<je; j++){
	for (i=is; i<ie; i++){
	  /* First coordinate */
	  U_plus = fmax(U[k][j][i],0.0); // max{U_{i-1/2,j,m},0.0} LHS boundary
	  U_minus = fmin(U[k][j][i+1],0.0); // min{U_{i+1/2,j,m},0.0} RHS boundary
#if !defined(SECOND_ORDER)
#ifdef SOLID_ANGLE
	  net_flux[k][j][i] += dt/(vol[k][j][i])*(dxa1[k]*dx2*(solid_U_plus[k]*I[k][j][i] + solid_U_minus[k]*I[k][j][i+1])-dxa1[k]*dx2*(solid_U_plus[k]*I[k][j][i-1] + solid_U_minus[k]*I[k][j][i]));
#else
#ifdef FACE_INTERPOLATE
	  //old solution to riemann problem: (i+1/2 net flux has + sign in divergence, i-1/2 has - sign)
	  /*	  if (fmin(0.0,U[k][j][i]) < 0) //left boundary wave travels left 
	    net_flux[k][j][i] -= dt/(vol[k][j][i])*(dxa1[k]*dx2*xFace_flux_l[k][j][i]);
	  if (fmax(0.0,U[k][j][i]) > 0) //left boundary wave travels right 
	    net_flux[k][j][i] -= dt/(vol[k][j][i])*(dxa1[k]*dx2*xFace_flux_r[k][j][i-1]);
	  if (fmin(0.0,U[k][j][i+1]) < 0) //right boundary wave travels left 
	    net_flux[k][j][i] += dt/(vol[k][j][i])*(dxa1[k]*dx2*xFace_flux_r[k][j][i+1]);
	  if (fmax(0.0,U[k][j][i+1]) > 0) //right boundary wave travels right 
	  net_flux[k][j][i] += dt/(vol[k][j][i])*(dxa1[k]*dx2*xFace_flux_l[k][j][i]); */

	  //divergence theorem
	  net_flux[k][j][i] += dt/(vol[k][j][i])*(dxa1[k]*dx2*xFace_flux[k][j][i+1] - dxa1[k]*dx2*xFace_flux[k][j][i]); 
#else
	  net_flux[k][j][i] += dt/(vol[k][j][i])*(dxa1[k]*dx2*(fmax(U[k][j][i+1],0.0)*I[k][j][i] + U_minus*I[k][j][i+1])-dxa1[k]*dx2*(U_plus*I[k][j][i-1] + fmin(U[k][j][i],0.0)*I[k][j][i]));
#endif //FACE_INTERPOLATE
#endif //SOLID_ANGLE
#else
	  /* Second order fluxes */
	  if (U[k][j][i+1] > 0.0){ //middle element is always the upwind element
	    imu[0] = I[k][j][i-1];
	    imu[1] = I[k][j][i];  
	    imu[2] = I[k][j][i+1];  
	    pr[0] = 0;
	    pr[1] = 0;
	    pr[2] = 0; 
	  }
	  else{
	    imu[0] = I[k][j][i+2];
	    imu[1] = I[k][j][i+1];  
	    imu[2] = I[k][j][i];  
	    pr[0] = 0;
	    pr[1] = 0;
	    pr[2] = 0; 
	  }
	  flux_limiter_x1_r= flux_PLM(dx1,imu);
	  flux_limiter_x1_r = flux_PLM_athena(pr, 2, dt, dx1, fabs(U[k][j][i+1]), imu);
#ifdef DEBUG
	  /* Debug: probe the limited fluxes */
	  if (n==n_debug && k==k_debug & j == j_debug && i == i_debug){
	    printf("-------------------\n");
	    printf("i+1/2 flux limiter \n");
	    printf("-------------------\n"); 
	    flux_limiter_x1_r = flux_PLM_athena_debug(pr, 2, dt, dx1, fabs(U[k][j][i+1]), imu);
	  }
#endif
	  //F^H_{i+1/2,j}
	  //net_flux[k][j][i] -= dt/(kappa[k][j][i]*dx1)*(x1_b[i+1]*(1-dt*fabs(U[k][j][i+1])/(dx1))*fabs(U[k][j][i+1])*flux_limiter_x1_r/2);

	  if (U[k][j][i] > 0.0){
	    imu[0] = I[k][j][i-2];  //points to the two preceeding bins; 
	    imu[1] = I[k][j][i-1];  
	    imu[2] = I[k][j][i];  
	    pr[0] = 0;
	    pr[1] = 0;
	    pr[2] = 0;
	  }
	  else{
	    imu[0] = I[k][j][i+1]; //centered around current bin
	    imu[1] = I[k][j][i];  
	    imu[2] = I[k][j][i-1];  
	    pr[0] = 0;
	    pr[1] = 0;
	    pr[2] = 0;
	  }
	  flux_limiter_x1_l= flux_PLM(dx1,imu);
	  flux_limiter_x1_l = flux_PLM_athena(pr, 2, dt, dx1, fabs(U[k][j][i]), imu);
#ifdef DEBUG
	  /* Debug: probe the limited fluxes */
      	  if (n==n_debug && k==k_debug & j == j_debug && i == i_debug){
	    printf("-------------------\n");
	    printf("i-1/2 flux limiter \n");
	    printf("-------------------\n"); 
	    flux_limiter_x1_l = flux_PLM_athena_debug(pr, 2, dt, dx1, fabs(U[k][j][i]), imu);
	    printf("limited x fluxes: i-1/2 = %.5e i+1/2 = %.5e \n",flux_limiter_x1_l,flux_limiter_x1_r);  
	    printf("dxa1[k]*dx2*flux_limiter_x1_r*U[k][j][i+1] = %.9e \n",dxa1[k]*dx2*flux_limiter_x1_r*U[k][j][i+1]);
	    printf("dxa1[k]*dx2*flux_limiter_x1_l*U[k][j][i] = %.9e \n",dxa1[k]*dx2*flux_limiter_x1_l*U[k][j][i]);
	    printf("dt/(vol[k][j][i]) = %.9e \n",dt/(vol[k][j][i])); 
	  }
#endif //DEBUG
	  
	  //F^H_{i-1/2,j}
	  //net_flux[k][j][i] += dt/(kappa[k][j][i]*dx1)*(x1_b[i]*(1-dt*fabs(U[k][j][i])/(dx1))*fabs(U[k][j][i])*flux_limiter_x1_l/2);
	  //	  net_flux[k][j][i] = dt/(kappa[k][j][i]*dx1)*(x1_b[i+1]*flux_limiter_x1_r*fabs(U[k][j][i+1]) - x1_b[i]*flux_limiter_x1_l*fabs(U[k][j][i]));

	  //working formula for copied flux_PLM
	  //	  net_flux[k][j][i] = dt/(kappa[k][j][i]*dx1)*(x1_b[i+1]*flux_limiter_x1_r*U[k][j][i+1] - x1_b[i]*flux_limiter_x1_l*U[k][j][i]);

	  //rewrite for bruls discretization
	  net_flux[k][j][i] += dt/(vol[k][j][i])*(dxa1[k]*dx2*flux_limiter_x1_r*U[k][j][i+1] - dxa1[k]*dx2*flux_limiter_x1_l*U[k][j][i]);

#endif
	  /* Second coordinate */
	  V_plus = fmax(V[k][j][i],0.0); // max{V_{i,j-1/2},0.0} LHS boundary
	  V_minus = fmin(V[k][j+1][i],0.0); // min{V_{i,j+1/2},0.0} RHS boundary

#if !defined(SECOND_ORDER)
#ifdef SOLID_ANGLE
	  net_flux[k][j][i] += dt/(vol[k][j][i])*(dx1*dxa1[k]*(solid_V_plus[k]*I[k][j][i] + solid_V_minus[k]*I[k][j+1][i])-dx1*dxa1[k]*(solid_V_plus[k]*I[k][j-1][i] + solid_V_minus[k]*I[k][j][i]));
#else
#ifdef FACE_INTERPOLATE
	  //old solution to riemann problem: (i+1/2 net flux has + sign in divergence, i-1/2 has - sign)
	  /*	  if (fmin(0.0,V[k][j][i]) < 0) //left boundary wave travels left 
	    net_flux[k][j][i] -= dt/(vol[k][j][i])*(dxa1[k]*dx1*yFace_flux_l[k][j][i]);
	  if (fmax(0.0,V[k][j][i]) > 0) //left boundary wave travels right 
	    net_flux[k][j][i] -= dt/(vol[k][j][i])*(dxa1[k]*dx1*yFace_flux_r[k][j-1][i]);
	  if (fmin(0.0,V[k][j+1][i]) < 0) //right boundary wave travels left 
	    net_flux[k][j][i] += dt/(vol[k][j][i])*(dxa1[k]*dx1*yFace_flux_r[k][j+1][i]);
	  if (fmax(0.0,V[k][j+1][i]) > 0) //right boundary wave travels right 
	  net_flux[k][j][i] += dt/(vol[k][j][i])*(dxa1[k]*dx1*yFace_flux_l[k][j][i]); */

	  //new application of divergence theorem
	  net_flux[k][j][i] += dt/(vol[k][j][i])*(dxa1[k]*dx1*yFace_flux[k][j+1][i] - dxa1[k]*dx1*yFace_flux[k][j][i]); 
#else
	  net_flux[k][j][i] += dt/(vol[k][j][i])*(dx1*dxa1[k]*(fmax(V[k][j+1][i],0.0)*I[k][j][i] + V_minus*I[k][j+1][i])-dx1*dxa1[k]*(V_plus*I[k][j-1][i] + fmin(V[k][j][i],0.0)*I[k][j][i]));
#endif //FACE_INTERPOLATE
#endif //SOLID_ANGLE 
#else
	  /* Second order fluxes */
	  if (V[k][j+1][i] > 0.0){
	    imu[0] = I[k][j-1][i];  //points to the two preceeding bins; 
	    imu[1] = I[k][j][i];  
	    imu[2] = I[k][j+1][i];  
	  }
	  else{
	    imu[0] = I[k][j+2][i];
	    imu[1] = I[k][j+1][i];  
	    imu[2] = I[k][j][i];  
	  }
	  flux_limiter_x2_r= flux_PLM(dx2,imu);
	  flux_limiter_x2_r = flux_PLM_athena(pr, 2, dt, dx2, fabs(V[k][j+1][i]), imu);//pr isnt dereferenced if dir!=1
#ifdef DEBUG	 
	  /* Debug: probe the limited fluxes */
	  if (n==n_debug && k==k_debug & j == j_debug && i == i_debug){
	    printf("-------------------\n");
	    printf("j+1/2 flux limiter \n"); 	   
	    printf("-------------------\n");
	    flux_limiter_x2_r = flux_PLM_athena_debug(pr, 2, dt, dx2, fabs(V[k][j+1][i]), imu); 
	  }
#endif 
	  //G^H_{i,j+1/2}
	  //	  net_flux[k][j][i] -= dt/(kappa[k][j][i]*dx2)*((1-dt*fabs(V[k][j+1][i])/(kappa[k][j][i]*dx2))*fabs(V[k][j+1][i])*flux_limiter_x2_r/2);
	  if (V[k][j][i] > 0.0){
	    imu[0] = I[k][j-2][i];  //points to the two preceeding bins; 
	    imu[1] = I[k][j-1][i];  
	    imu[2] = I[k][j][i];  
	  }
	  else{
	    imu[0] = I[k][j+1][i]; //centered around current bin
	    imu[1] = I[k][j][i];  
	    imu[2] = I[k][j-1][i];  
	  }
	  flux_limiter_x2_l = flux_PLM(dx2,imu);
	  flux_limiter_x2_l = flux_PLM_athena(pr, 2, dt, dx2, fabs(V[k][j][i]), imu);
#ifdef DEBUG
	  /* Debug: probe limited slopes */
	  if (n==n_debug && k==k_debug & j == j_debug && i == i_debug){
	    printf("-------------------\n");
	    printf("j-1/2 flux limiter \n");
	    printf("-------------------\n"); 	   
	    flux_limiter_x2_l = flux_PLM_athena_debug(pr, 2, dt, dx2, fabs(V[k][j][i]), imu); 
	    printf("limited y fluxes: j-1/2 = %.5e j+1/2 = %.5e \n",flux_limiter_x2_l,flux_limiter_x2_r); 
	  }
#endif

	  //G^H_{i,j-1/2}
	  //	  net_flux[k][j][i] += dt/(kappa[k][j][i]*dx2)*((1-dt*fabs(V[k][j][i])/(kappa[k][j][i]*dx2))*fabs(V[k][j][i])*flux_limiter_x2_l/2);

	  //working formula for copied flux_PLM
	  //	  net_flux[k][j][i] += dt/(kappa[k][j][i]*dx2)*(flux_limiter_x2_r*V[k][j+1][i] - flux_limiter_x2_l*V[k][j][i]);

	  //rewrite for bruls discretization
	  net_flux[k][j][i] += dt/(vol[k][j][i])*(dx1*dxa1[k]*flux_limiter_x2_r*V[k][j+1][i] - dxa1[k]*dx1*flux_limiter_x2_l*V[k][j][i]);
#endif
	}
      }
    }	
    
    /* Apply fluxes */
    for (k=ks; k<ke; k++){
      for (j=js; j<je; j++){
	for (i=is; i<ie; i++){
	  I[k][j][i] -= net_flux[k][j][i];	
	}
      }
    }
    
    /* Add source terms */
      for (k=ks; k<ke; k++){
	for (j=js; j<je; j++){
	  for (i=is; i<ie; i++){
	    I[k][j][i] += source[k][j][i]*I[k][j][i]*dt;
	  }
	}
      }

      /* Reset any fixed points inside domain for outputting purposes */
    for (k=ks; k<ke; k++){
      for (j=js; j<je; j++){
	for (i=is; i<ie; i++){
#ifdef TEST5
	  if (i == is + nx1_r/2 && j == js + nx2_r/2 && k >= source_k_start  && k <= source_k_end){ //integer division
	    I[k][j][i] = 1.0;	   
	  }
#endif 
	  fixed_I = bc_interior(x1[i], x2[j], xa1[k], n*dt);
	  if (fixed_I != 0.0)
	    I[k][j][i] = fixed_I; 
	}
      }
    }


    /* Output */
    //for now, explicitly copy subarray corresponding to real zonal inf (as opposed to setting pointers):
    index=0; 
    indexJ=0;
    for (j=0; j<nx2_r; j++){
      for (i=0; i<nx1_r; i++){
	realJ[indexJ] = 0;
	indexJ++;
      }
    }

    for (k=ks; k<ke; k++){
      indexJ =0; 
      for (j=js; j<je; j++){
	for (i=is; i<ie; i++){
#ifdef DEBUG
	  /* Debug: probe the intensity of a cell */
	  if (n==n_debug && k==k_debug & j == j_debug && i == i_debug){
	    printf("(%d,%d,%d) (U,V) = (%lf, %lf) \n",k,j,i,U[k][j][i],V[k][j][i]);
	    printf("Current I = %.9e Previous I = %.9e Net Flux = %.9e\n",I[k][j][i],I[k][j][i] +net_flux[k][j][i], net_flux[k][j][i]); 
	  }
#endif	
	  //index =(j-num_ghost)*nx2_r + (i-num_ghost); 
	  realI[index] = (float) I[k][j][i];
	  /* Diagonal beam test only: check for nonphysical extrema on x and y edges of beam */
	  if ((k==1 && i > source_i && j== source_j && i<ie))// || (j > source_j && i== source_i && j<je))
	    if (I[k][j][i+1] > I[k][j][i]){
	      printf("Nonphysical x-extrema detected! n,k,j,i = (%d,%d,%d,%d) \n",n,k,j,i); 
	      printf("Current I = %.5e Previous I = %.5e Net Flux = %.5e\n",I[k][j][i],I[k][j][i] +net_flux[k][j][i], net_flux[k][j][i]); 
	      printf("Spatially adjacent intensities (prev step):\n I_{i-2} = %.5e I_{i-1} = %.5e I_{i+1} = %.5e I_{i+2} = %.5e \n I_{j-2} = %.5e I_{j-1} = %.5e I_{j+1} = %.5e I_{j+2} = %.5e \n",I[k][j][i-2] +net_flux[k][j][i-2],I[k][j][i-1] +net_flux[k][j][i-1],I[k][j][i+1] +net_flux[k][j][i+1],I[k][j][i+2] +net_flux[k][j][i+2],I[k][j-2][i] +net_flux[k][j-2][i],I[k][j-1][i] +net_flux[k][j-1][i],I[k][j+1][i] +net_flux[k][j+1][i],I[k][j+2][i] +net_flux[k][j+2][i]); 
	      printf("Spatially adjacent intensities (current step):\n I_{i-2} = %.5e I_{i-1} = %.5e I_{i+1} = %.5e I_{i+2} = %.5e \n I_{j-2} = %.5e I_{j-1} = %.5e I_{j+1} = %.5e I_{j+2} = %.5e \n",I[k][j][i-2],I[k][j][i-1],I[k][j][i+1],I[k][j][i+2],I[k][j-2][i],I[k][j-1][i],I[k][j+1][i],I[k][j+2][i]);
	    }  
	  /*trap any negative intensity */
	  if (I[k][j][i] <0.0){
	    printf("Negative intensity found at  n,k,j,i = (%d,%d,%d,%d) \n",n,k,j,i); 
	    printf("Current I = %.5e Previous I = %.5e Net Flux = %.5e\n",I[k][j][i],I[k][j][i] +net_flux[k][j][i], net_flux[k][j][i]); 
	    printf("Edge velocities: U_{i-1/2} = %lf U_{i+1/2} = %lf V_{j-1/2} = %lf V_{j+1/2} = %lf \n",U[k][j][i],U[k][j][i+1],V[k][j][i],V[k][j][i+1]); 
	    printf("Spatially adjacent intensities (prev step):\n I_{i-2} = %.5e I_{i-1} = %.5e I_{i+1} = %.5e I_{i+2} = %.5e \n I_{j-2} = %.5e I_{j-1} = %.5e I_{j+1} = %.5e I_{j+2} = %.5e \n",I[k][j][i-2] +net_flux[k][j][i-2],I[k][j][i-1] +net_flux[k][j][i-1],I[k][j][i+1] +net_flux[k][j][i+1],I[k][j][i+2] +net_flux[k][j][i+2],I[k][j-2][i] +net_flux[k][j-2][i],I[k][j-1][i] +net_flux[k][j-1][i],I[k][j+1][i] +net_flux[k][j+1][i],I[k][j+2][i] +net_flux[k][j+2][i]); 
	    printf("Spatially adjacent intensities (current step):\n I_{i-2} = %.5e I_{i-1} = %.5e I_{i+1} = %.5e I_{i+2} = %.5e \n I_{j-2} = %.5e I_{j-1} = %.5e I_{j+1} = %.5e I_{j+2} = %.5e \n",I[k][j][i-2],I[k][j][i-1],I[k][j][i+1],I[k][j][i+2],I[k][j-2][i],I[k][j-1][i],I[k][j+1][i],I[k][j+2][i]);
	    return(1); 
	  }
	  /*compute zeroth angular moment of the radiation field*/
	  realJ[indexJ] += (float) I[k][j][i]*pw[k];
	  indexJ++; 
	  index++;
	}
      }
    }
    /* Copy the radiation energy field from first angular bin to all angular bins */
    indexJ=nx1_r*nx2_r; 
    for (k=ks+1; k<ke; k++){
      for (j=0; j<nx2_r; j++){
	for (i=0; i<nx1_r; i++){
	  realJ[indexJ] = realJ[j*nx1_r+i]; 
	  indexJ++;
	}
      }
    }
   
    sprintf(filename,"rte-%.3d.vtk",n); 
    if(!OUTPUT_INTERVAL){
      if (n==nsteps-1) //for only the final result
	write_curvilinear_mesh(filename,3,dims, pts, nvars,vardims, centering, varnames, vars);}
    else{
      if (!(n%OUTPUT_INTERVAL)) 
	write_curvilinear_mesh(filename,3,dims, pts, nvars,vardims, centering, varnames, vars);}
      
      printf("step: %d time: %lf max{I} = %0.7lf min{I} = %0.7lf sum{I} = %0.7lf \n",
	     n,n*dt,find_max(realI,nxa1*nx1_r*nx2_r),find_min(realI,nxa1*nx1_r*nx2_r),sum(realI,nxa1*nx1_r*nx2_r));
  }

/* Compute J error from analytic solution */
#ifdef ANALYTIC_SOLUTION
  double error_L1, error_L2, error_max; 
  error_L1 = 0.0; 
  error_L2 = 0.0; 
  error_max = 0.0;
  temp = 0.0;  
  indexJ = 0;
  //only use points inside first spatial quadrant-- not anymore! 
  for (j=0; j<nx2_r; j++){ //even though realJ is a flattened 3D structure, we can reference only the first z height
    for (i=0; i<nx1_r; i++){ 
      if (!(i+num_ghost == source_i && j+num_ghost == source_j)){
	temp = fabs(realJ[indexJ] -realJ_analytic[indexJ]);
	error_L1 += temp;
	error_L2 += pow(realJ[indexJ] -realJ_analytic[indexJ],2);
	if (temp > error_max)
	  error_max = temp; 
      }
      /*	else { //only turn on for diagonal test
		if (realJ[indexJ] != 0.0 || realJ_analytic[indexJ] != 0.0){
		printf("nonzero J outside of first quadrant! \n"); 
		printf("%d, %d numericalJ = %.5e analyticJ = %.5e \n",i,j,realJ[indexJ],realJ_analytic[indexJ]); }
		} */
      indexJ++;
    }
  }
  error_L2 = sqrt(error_L2); 
  printf("Errors (L1,L2,MAX) = (%lf, %lf, %lf) \n",error_L1,error_L2,error_max);  
#endif// ANALYTIC_SOLUTION
  return(0); 
}

/* Map to physical (cartesian) coordinates */
double X_physical(double x1, double x2){
  double x_cartesian = x1*cos(x2); 
  return(x_cartesian); 
}
double Y_physical(double x1, double x2){
  double y_cartesian = x1*sin(x2); 
  return(y_cartesian); 
}

double initial_condition(double x, double y, double xa1){
  return(0.0);
}

double bc_x1i(double x, double y, double xa1, double t){
  return(0.0);
}
//bc at outermost radius
double bc_x1f(double x, double y, double xa1, double t){
  return(0.0);
}
//this bc is specified in terms of polar coordinates
float bc_x1f_polar(double phi, double xa1, double t){
  return(0.0);
}
//if you want to fix a point inside the domain at a particular intensity 
double bc_interior(double r, double phi, double xa1, double t){
  return(0.0);
}
//bc at phi=0.0
double bc_x2i(double x, double y, double xa1, double t){
  return(0.0);
}
//bc at phi_final
double bc_x2f(double x, double y, double xa1, double t){
  return(0.0);
}

float find_max(float a[], int n) {
  int i,index;
  float max; 
  max = a[0];
  index = 0;
  for (i = 1; i < n; i++) {
    if (a[i] > max) {
      index = i;
      max = a[i];
    }
  }
  return(max); 
}

double find_max_double(double a[], int n) {
  int i,index;
  double max; 
  max = a[0];
  index = 0;
  for (i = 1; i < n; i++) {
    if (a[i] > max) {
      index = i;
      max = a[i];
    }
  }
  return(max); 
}

float find_min(float a[], int n) {
  int i,index;
  float min; 
  min = a[0];
  index = 0;
  for (i = 1; i < n; i++) {
    if (a[i] < min) {
      index = i;
      min = a[i];
    }
  }
  return(min); 
}

float sum(float a[], int n) {
  int i,index;
  float sum; 
  sum = a[0];
  for (i = 1; i < n; i++) {
    sum+= a[i]; 
  }
  return(sum); 
}


int uniform_angles2D(int N, double phi_z, double *pw, double **mu, double **mu_b, double *xa1,double *xa1_b){
/* Generate uniform 2D discretization of mu */

  //   Input: N-- Number of mu level cells in one \hat{k}^i dimension
  //  (note, there are n+1, n/2+1 boundaries, including -1,1 and 0,1) N MUST BE EVEN, should be greater than 6
  //   phi: angle down from the z-axis that fixes the plane of propagation directions. pi/2 for k_z=0

  //   Output: 
  //   nxa: angular parameter coordinate dimensions. currently, it is 2-dim with nxa(1) = N nxa(2) = N +1.   

  // The algorithm discretizes theta \in [0, 2pi ) then generates the boundary rays mu_b from these \theta
  // The actual sampling rays are generated by averaging the \theta of the boundary rays and then 
  // using spherical polar coordinates to convert to k^x, k^y, k^z

  // Keep data structures same as 3D case so that 2D problems are a subset									  
  int i,j,k;
  int nxa1 = N;
  double  dxa1 = 2*M_PI/(nxa1); //dont mesh all the way to 2pi
  xa1_b[0] = 0.0;
  xa1[0] = dxa1/2;
  for (i=1; i<nxa1; i++){
    xa1_b[i] = xa1_b[i-1] + dxa1; 
    xa1[i] = xa1[i-1] + dxa1; 
  }

  //add another boundary ray for purposes of VTK output
  xa1_b[nxa1]  = xa1_b[nxa1-1] + dxa1; 
 
  for (i=0; i<= nxa1; i++){
    //  for (i=0; i< nxa1; i++){
    mu_b[i][0] = cos(xa1_b[i])*sin(phi_z);
    mu_b[i][1] = sin(xa1_b[i])*sin(phi_z);
    mu_b[i][2] = cos(phi_z);
  }
  //periodicity of the domain implies duplication of 1 theta bin first and last column/row are identical for kx, ky
  for (i=0; i<nxa1; i++){
    mu[i][0] = cos(xa1[i])*sin(phi_z);
    mu[i][1] = sin(xa1[i])*sin(phi_z);
    mu[i][2] = cos(phi_z);
  }

  /*------------------CALCULATE QUADRATURE WEIGHTS ------------------- */
  // Should this be proportional to the "size of the cell"?
  for (i=0; i<nxa1; i++)
    pw[i] = 1.0/(nxa1);
  return(nxa1);  
}

/* Function for returning pointer that allows A[i][j] access of array 
   such that i increases first in memory, then second index */
// pass in contiguous chunk of memory of size n1*n2*n3
// creates many excess pointers ....
double **allocate_2D_contiguous(double *a, int n1, int n2){ //this is reordered since advection
  int i; 
  double **a_p = (double **) malloc(sizeof(double *)*n1);
  for(i=0; i<n1; i++)
    a_p[i] = malloc(sizeof(double)*n2); //&(a[n2*i]);    
  return(a_p); 
}

/* Function for returning pointer that allows A[][][] access of array */
double ***allocate_3D_contiguous(double *a, int n1, int n2, int n3){
  int i,j; 
  double ***a_p = (double ***) malloc(sizeof(double **)*n1);
  double **a_p2 = (double **) malloc(sizeof(double *)*n2*n1);
  for(i=0; i<n1; i++){
    a_p[i] = (double **) malloc(sizeof(double *)*n2); //&(a_p2[n2*n3*i]);
    for (j=0; j< n2; j++){
      a_p[i][j] = malloc(sizeof(double)*n3); //&(a[n2*n3*i + n3*j]);    
    }
  } 
  return(a_p);
}
/* a near duplication of the function in ATHENA FullRT_flux.c */
double flux_PLM(double ds,double *imu){
  //the upwind slope
  double delq1 = (imu[2] - imu[1])/ds;
  double delq2 = (imu[1] - imu[0])/ds;
  double dqi0;
  if(delq1*delq2 >0.0) //minmod function
    dqi0 = 2.0*delq1*delq2/(delq1+delq2);
  else
    dqi0 = 0.0;
  //unknown why ds is in this function. might have to do with nonuniform grids
  return(ds*dqi0); 
}

double flux_PLM_athena(double r[3], int dir, double dt, double ds, double vel, double imu[3])
{
  /* for each ray, we always use the upwind specific intensity , therefore velocity is always positive */
  /* imup[0:2] is i-2, i-1, i*/

  double dqi0, delq1, delq2;
  double *imup;
  double distance;
  double *pr;
  pr = &(r[2]);
  double geom1 = 1.0, geom2 = 1.0;
  double geom3 = 1.0, geom4 = 1.0; /* geometric weighting factor */

  imup = &(imu[2]);

  /* The upwind slope */
  delq1 = (imup[0] - imup[-1]) / ds;
  delq2 = (imup[-1] - imup[-2]) / ds;
  /* Only need to apply the weighting factor when we want to calculate
   * flux along radial direction */
  /*  if (dir ==1 ){
    geom1 = 1.0/(1.0 - ds * ds/(12.0 * pr[0] * pr[-1]));
    geom2 = 1.0/(1.0 - ds * ds/(12.0 * pr[-1] * pr[-2]));
    
    delq1 *= geom1;
    delq2 *= geom2;
    } */

  if(delq1 * delq2 > 0.0){
    /* Original ATHENA source code that may result in denormal numbers without flooring */
    //    dqi0 = 2.0 * delq1 * delq2 / (delq1 + delq2);
    dqi0 = 2.0 * delq2 / (delq1 + delq2);
    dqi0 *= delq1; 
  }
  else{
    dqi0 = 0.0;
  }

  /* Take into account the curvature effect for time averaged interface state */
  /* See eq. 64 of skinner & ostriker 2010 */
  /*  if (dir == 1){
    geom3 = 1.0 - ds /(6.0 * pr[-1]);
    geom4 = 1.0 - vel * dt/(6.0*(0.5 * (pr[-1] + pr[0]) - 0.5 * vel * dt));
    } */

  distance = ds * geom3;
  distance -= ((vel * dt) * geom4);

  /* The upwind flux */
  return(imup[-1] + distance * (dqi0/2.0)); 
}

double flux_PLM_athena_debug(double r[3], int dir, double dt, double ds, double vel, double imu[3])
{ //modified for printing out every step of process
  /* for each ray, we always use the upwind specific intensity , therefore velocity is always positive */
  /* imup[0:2] is i-2, i-1, i*/

  double dqi0, delq1, delq2;
  double *imup;
  double distance;
  double *pr;
  pr = &(r[2]);
  double geom1 = 1.0, geom2 = 1.0;
  double geom3 = 1.0, geom4 = 1.0; /* geometric weighting factor */

  imup = &(imu[2]);

  /* The upwind slope */
  delq1 = (imup[0] - imup[-1]) / ds;
  delq2 = (imup[-1] - imup[-2]) / ds;
  /* Only need to apply the weighting factor when we want to calculate
   * flux along radial direction */
  /*  if (dir ==1 ){
    geom1 = 1.0/(1.0 - ds * ds/(12.0 * pr[0] * pr[-1]));
    geom2 = 1.0/(1.0 - ds * ds/(12.0 * pr[-1] * pr[-2]));
    
    delq1 *= geom1;
    delq2 *= geom2;
    } */

  if(delq1 * delq2 > 0.0){
    /* Original ATHENA source code that may result in denormal numbers without flooring */
    //    dqi0 = 2.0 * delq1 * delq2 / (delq1 + delq2);
    dqi0 = 2.0 * delq2 / (delq1 + delq2);
    dqi0 *= delq1;  
  }
  else{
    dqi0 = 0.0;
  }

  /* Take into account the curvature effect for time averaged interface state */
  /* See eq. 64 of skinner & ostriker 2010 */
  /*  if (dir == 1){
    geom3 = 1.0 - ds /(6.0 * pr[-1]);
    geom4 = 1.0 - vel * dt/(6.0*(0.5 * (pr[-1] + pr[0]) - 0.5 * vel * dt));
    } */

  distance = ds * geom3;
  distance -= ((vel * dt) * geom4);
  //  printf("(imup[0] - imup[-1]) = %.5e \n",(imup[0] - imup[-1]));
  printf("delq1 * delq2 = %.5e \n",delq1 * delq2);
  printf("2.0 * delq1 * delq2 = %.5e \n",2.0 * delq1 * delq2);
  printf("delq1 + delq2 = %.5e\n",delq1 + delq2); 
  printf("dqi0*ds = %.5e \n",dqi0*ds);
  printf("delq1 = %.5e delq2 = %.5e dqi0 = %.5e upwind flux = %.5e\n",delq1,delq2,dqi0,imup[-1] + distance * (dqi0/2.0)); 
  printf("ds = %.5e vel*dt = %.5e \n",ds,vel*dt); 
  printf("distance = %.5e, distance*dqi0/2.0 = %.5e imup[-1] = %.5e\n",distance,distance*dqi0/2.0,imup[-1]); 
  /* The upwind flux */
  return(imup[-1] + distance * (dqi0/2.0)); 
  // the issue in the current 128x128x196 diagonal test case is that the PL reconstruction is producing a value that is negative at I_{j-1/2}
  // at (61,0,101,92). When it should be between 2.45658e-167 and 0. The monotonicity requirement on the limited slope dqi0 should prevent this overshoot

  /* delq1 = -3.11986e-165 delq2 = -1.20601e-159 dqi0 = -8.19339e-165 upwind flux = -7.54539e-168
distance = 7.83832e-03, distance*dqi0/2.0 = -3.21112e-167 imup[-1] = 2.45658e-167
y flux limiter: j-1/2 = -7.54539e-168 j+1/2 = 0.00000e+00  */

  /* Current I = -2.51040e-170 Previous I = 0.00000e+00 Net Flux = 2.51040e-170
Edge velocities: U_{i-1/2} = 0.999829 U_{i+1/2} = 0.999829 V_{j-1/2} = 0.016027 V_{j+1/2} = 0.016027 
Spatially adjacent intensities (prev step):
 I_{i-2} = 3.99954e-168 I_{i-1} = 3.21844e-170 I_{i+1} = 0.00000e+00 I_{i+2} = 0.00000e+00 
 I_{j-2} = 9.49614e-162 I_{j-1} = 2.45658e-167 I_{j+1} = 0.00000e+00 I_{j+2} = 0.00000e+00 
Spatially adjacent intensities (current step):
 I_{i-2} = 6.18134e-167 I_{i-1} = 1.19385e-168 I_{i+1} = 0.00000e+00 I_{i+2} = 0.00000e+00 
 I_{j-2} = 4.15855e-160 I_{j-1} = 3.89846e-166 I_{j+1} = 0.00000e+00 I_{j+2} = 0.00000e+00 */

  //turning CFL down to 0.3 results in a slightly earlier negative intensity
  //  (46,0,103,69) 

  //turning CFL up to 0.8 results in many nonphysical extrema
  // and much earlier neagitve intensity (24,19,66,82) 

  // the velocity should not be an issue since it is constant and positive everywhere. 
  // the y-edge velocity is small 0.016027, which makes it hard for me to believe that this is a cfl violation
}
int bruls_angles2D(int N, double phi_z, double *pw, double **mu, double **mu_b, double *xa1,double *xa1_b){
  int num_rays = N*(N+2)/2; //up to 84 in 2D
  int num_rays_per_quadrant = num_rays/4; 


  //made modifications relative to MATLAB code to force phi_z=pi
  double *ordinates = malloc(sizeof(double)*N/2); 
  ordinates[0] = sqrt(1.0/(3.0*(N-1))); //changed
  double dcos = (1-3*pow(ordinates[0],2))*2/(N-2); //changed 

  int i, j,k,l,m;
  /* Choice of first cosine is arbitrary. In analogy to Gaussian quadrature */
  for (i=0; i<N/2; i++){
    ordinates[i] = sqrt(ordinates[0]*ordinates[0] + i*dcos); 
  }

  // Derive weights for integration. Ref: Bruls 1999 appendix
  double *W = malloc(sizeof(double)*N/2); 
  for (i=0; i<N/2; i++){
    W[i] = sqrt(4.0/(3.0*(N-1)) + i*2.0/(N-1));  //should I change?
    //    printf("W[%d] = %lf\n",i,W[i]); 
  }
  //compute level weights
  double *lw = malloc(sizeof(double)*N/2); 
  lw[0] = W[0];
  double sum =lw[0]; 
  for (i=1; i<N/2-1; i++){
    lw[i] = W[i] - W[i-1];
    sum += lw[i]; 
  }

  /* following ATHENA, we correct the last level weight to renormalize...not
     sure why this happens */
  lw[N/2-1] = 1.0 - sum; 
  
  /* Permutation matrix */
  double **pmat; 
  double *data_pmat = malloc(sizeof(double)*N/2*N/2); 

  pmat = allocate_2D_contiguous(data_pmat,N/2,N/2); 
  int *plab; 
  plab = malloc(sizeof(int *)*num_rays); 
  
  int **pl = malloc(sizeof(int *)*N/2); 
  //initialize pl to zeros
  for (i=0; i<N/2; i++)
    pl[i] = malloc(sizeof(int)*3); 
  /*  for (i=0; i<N/2; i++)
    for (j=0; j<3; j++)
      pl[i][j] = 0;  */

  int ray_count = 0; 
  int np =0; 
  int ip =0; 
  /* To select valid rays for the quadrature, we follow the index selection
  rule: sigma = i + j + k
  i.e. indices of the direction cosines must equal ntheata/2 +2
  for proper normalization in 3D space */
  for (i=0; i<N/2; i++){
    for (j=0; j<N/2-i; j++){
      k = N/2-1-i-j; //assume this is correct 
      mu[ray_count][0] = ordinates[i]; 
      mu[ray_count][1] = ordinates[j]; 
      ip = permutation(i,j,k,pl,np); 
      /*      printf("ip = %d np = %d i,j,k = %d,%d,%d\n",ip,np, i,j,k); 
      printf("pl =\n");
      for (l=0; l<N/2; l++){
	for (m=0; m<3; m++)
	  printf("%d ", pl[l][m]); 
	printf("\n");
	} */
      if (ip == -1){ //ray indices havent been loaded yet
	pl[np][0] = i; 
	pl[np][1] = j; 
	pl[np][2] = k; 
	pmat[i][np]++; 
	plab[ray_count] = np; 
	np = np+1; 
      }
      else {
	pmat[i][ip]++;
	plab[ray_count] = ip; 	     
      }
      ray_count++;  
    }
  }
  assert(ray_count == num_rays_per_quadrant); 

  /* discretization symmetry: reflect across second axis */
  for (i=0; i<ray_count; i++){
    mu[ray_count+i][0] = -mu[i][0]; 
    mu[ray_count+i][1] = mu[i][1]; 
  }
  ray_count*=2;  
  /* discretization symmetry: reflect across second axis */
  //flip the order for this so that the rays are in order of xa1
  for (i=0; i<ray_count; i++){
    mu[2*ray_count-(i+1)][0] = mu[i][0]; 
    mu[2*ray_count-(i+1)][1] = -mu[i][1]; 
  }
  ray_count*=2; 

  /* Solve system of equations to calculate families of point weights */
  double *wpf = malloc(sizeof(double)*(N/2-1)); 
  /*  printf("pmat=\n");
  for (l=0; l<N/2; l++){
    for (m=0; m<N/2; m++)
      printf("%lf ", pmat[l][m]); 
    printf("\n");
  }
  printf("lw=\n");
  for (l=0; l<N/2; l++)
  printf("%lf ", lw[l]);  */

  gaussianelim(pmat, lw, wpf, N/2-1, 1);
  /*  for (i=0; i<N/2-1; i++)
    printf("wpf[%d] = %lf\n",i,wpf[i]);
  for (i=0; i<ray_count; i++)
    printf("plab[%d] = %d\n",i,plab[i]);
  */
  for(i=0; i<num_rays_per_quadrant; i++){
    pw[i] = wpf[plab[i]]/4; 
    pw[i+num_rays_per_quadrant] = wpf[plab[i]]/4; 
    pw[i+2*num_rays_per_quadrant] = wpf[plab[i]]/4; 
    pw[i+3*num_rays_per_quadrant] = wpf[plab[i]]/4; 
  }
  
  /* Compute angles corresponding to sample rays */
  for (i=0; i<num_rays; i++){
    xa1[i] = atan2(mu[i][1],mu[i][0]) + M_PI;
  }
  /* Sort by increasing xa1 */
  qsort(xa1, num_rays, sizeof(double), cmpfunc); 
  //better way of resorting mu?
  //recompute mu
  for (i=0; i < num_rays; i++){
    mu[i][0] = cos(xa1[i])*sin(phi_z); 
    mu[i][1] = sin(xa1[i])*sin(phi_z);
    mu[i][2] = cos(phi_z); 
  }

  /* Boundary angles lie halfway between adjacent sample rays */
  //i am assuming they are ordered in xi here, and:
  // the first ray is in the first quadrant
  // the last ray is in the last quadrant 
  xa1_b[0] = fmod((xa1[0] +2*M_PI + xa1[num_rays-1])/2, 2*M_PI); 
  for (i=1; i<num_rays; i++){
    xa1_b[i] = (xa1[i] + xa1[i-1])/2; 
  }
  //  xa1_b[num_rays] = xa1_b[0];  
  xa1_b[num_rays] = 2*M_PI; //can i assume this?
  for (i=0; i <= num_rays; i++){
    mu_b[i][0] = cos(xa1_b[i])*sin(phi_z); 
    mu_b[i][1] = sin(xa1_b[i])*sin(phi_z);
    mu_b[i][2] = cos(phi_z); 
  }
  return(num_rays); 
}
//stolen from ATHENA to calculate point weights
int permutation(int i, int j, int k, int ** pl, int np){
  int ip = -1;
  int l,m,n,o; 
  for (l=0; l<np; l++){
    for(m=0; m<3; m++){
      if (i == pl[l][m])
	for (n=0; n<3; n++){
	  if (n != m){
	    if (j == pl[l][n]){
	      for (o=0; o<3; o++){
		if(o!= m && o!= n){
		  if (k == pl[l][o]){
		    ip = l; 
		  }
		}
	      }
	    }
	  }
	}
    }
  }           
  return(ip); 
}

void gaussianelim(double **A, double *b, double *x, int n, int pivot) { //pivot 1 for partial row, 0 for no pivoting
  int i, j, pivot_pos, column;
  double **Ab; 
  double ratio,tmp;

  /* Allocate the new augmented matrix */
  Ab = (double **) malloc(sizeof(double)*(n));//allocate n pointers to rows of pointers
  for (i=0; i<n; i++){
    Ab[i] = (double *) malloc(sizeof(double)*(n+1));
  }

  /*Copy the old matrix and rhs into augmented matrix */
  for (i=0; i <n; ++i){     //rows
    for (j=0; j <n; ++j){ //collumns
      Ab[i][j] = A[i][j];
    }
  }
  for (i=0; i<n; i++) Ab[i][n] = b[i];

  /*Gaussian Elimination */
  /*Pivoting */
  double *temp;
  temp = (double *) calloc(sizeof(double),n+1);
  if (pivot) {
    for (j=0; j<n; j++){
      for (i=j+1; i<n; i++){ //check all rows beneath current pivot pos
	if (fabs(Ab[i][j]) > fabs(Ab[j][j])) {
	  temp = Ab[i];
	  Ab[i] = Ab[j];
	  Ab[j] = temp;
	}
      }
    }
  }
  
  /* Elimination of variables via basic row operations */
  for (i =0; i < (n-1); i++) {           //starting with first row, first column
    for (j=i+1; j<n; ++j) {              // take next row, first column
      ratio = Ab[j][i] / Ab[i][i];       //compute their ratio
      for (column = i; column <n+1; column++){    // and eliminate it from the second row
	Ab[j][column] -= (ratio * Ab[i][column]); // by subtracting it from all coef
      } //do this for all rows until the last row
    }
  }
  
  /* Back Subsitution */
  for (i=n-1; i>=0; i--){ //rows, starting from bottom
    tmp = 0; //Build a temporary dbl precision that holds the solved variables* the rows coefficients to subtract from RHS for each line
    for (j= i+1; j<n; j++) { //this inner loop doesnt do anything for the bottom row, since the equation req no more info
      tmp += x[j]*Ab[i][j];
    }
    x[i] = (Ab[i][n] - tmp)/ Ab[i][i]; //RHS - tmp (known variables, scaled) / unknown coeff
  }
}

int cmpfunc (const void * a, const void * b)
{
  double a_r = *(double*)a;
  double b_r = *(double*)b;
  if (a_r > b_r ){
    return(1);
  }
  else if (b_r > a_r){
    return(-1);
  }
  else {
    return(0);
  }
}

double min_cosine (double x, void *params){
  return(fmin(cos(x),0.0)); 
}
double max_cosine (double x, void *params){
  return(fmax(cos(x),0.0)); 
}
double min_sine (double x, void *params){
  return(fmin(sin(x),0.0)); 
}
double max_sine (double x, void *params){
  return(fmax(sin(x),0.0)); 
}
