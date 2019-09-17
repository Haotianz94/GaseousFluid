#ifndef _STDAFXNUMGRIDH_
#define _STDAFXNUMGRIDH_

#include <iostream>


#define CONNECTED 1
#define OBSTACLE 1
//#define OUTPUT 1
#define GAUSS_SEIDEL 1


#define NUMGRIDW 200
#define NUMGRIDH 200
#define CubeLength 1.0
//#define VISBLEW 1200
#define GRIDSIZE 6
#define DIFFUSION 0.001
#define VISCOSITY 0.001
#define TIMESTEP 0.01
#define ITERATION 30
#define FRAMERATE 32
#define DRAGSCALE 100
#define FLOWTIME 1
#define DENSITY 300
#define SPEED 4000
#define OBSTACLEX 30
#define LICL 30
//The viscosity matters a lot, when you decrease visc, you should also decrease speed


#ifdef CONNECTED
	#define IX(x, y) ( (x) == 0? NUMGRIDW + (y) * (NUMGRIDW+2) : ((x) == NUMGRIDW+1? 1 + (y) * (NUMGRIDW+2) : (x) + (y) * (NUMGRIDW+2)) )
#else	
	#define IX(x, y) ( (x) + (y) * (NUMGRIDW+2) )
#endif
#define IX2(x, y) ( (x) + (y) * (NUMGRIDW*GRIDSIZE+2) )
#define BOUNDED(x, y) ( ((x) < 1 || (x) > NUMGRIDW+1 || (y) < 1 || (y) > NUMGRIDH+1)? false : true)
	
#define PI 3.14159265
#define LENGTH _N*GRIDSIZE
#define SWAP(x0, x) {float *tmp = x0; x0 = x; x = tmp;}

//system output defines
#define REPORT(X) std::cout << #X << ": " << X << std::endl
#define PRINT(X) std::cout << X << std::endl
#define PRINT_ONELINE(X) std::cout << X

#endif

//Best for vortex street connected
/*
	#define NUMGRIDW 600
	#define NUMGRIDH 100
	#define CubeLength 1.0
	#define VISBLEW 1200
	#define GRIDSIZE 3
	#define DIFFUSION 0.01
	#define VISCOSITY 0.01
	#define TIMESTEP 0.01
	#define ITERATION 30
	#define FRAMERATE 32
	#define DRAGSCALE 100
	#define FLOWTIME 10
	#define DENSITY 100
	#define SPEED 20000
	#define OBSTACLEX 30
	#define LICL 30
*/

//Best for vortex street not connected
/*
	#define NUMGRIDW 150
	#define VISBLEW 1500
	#define NUMGRIDH 30
	#define GRIDSIZE 10
	#define DIFFUSION 0.01
	#define VISCOSITY 0.01
	#define TIMESTEP 0.01
	#define ITERATION 10
	#define FRAMERATE 32
	#define DRAGSCALE 100
	#define FLOWTIME 20
	#define DENSITY 100
	#define SPEED 200
*/

//Best for gas
/*
	#define NUMGRIDW 200
	#define NUMGRIDH 200
	#define CubeLength 1.0
	//#define VISBLEW 1200
	#define GRIDSIZE 3
	#define DIFFUSION 0.001
	#define VISCOSITY 0.001
	#define TIMESTEP 0.01
	#define ITERATION 30
	#define FRAMERATE 32
	#define DRAGSCALE 100
	#define FLOWTIME 1
	#define DENSITY 300
	#define SPEED 4000
	#define OBSTACLEX 30
	#define LICL 30
*/