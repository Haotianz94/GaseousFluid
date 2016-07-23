#ifndef _STDAFX_H_
#define _STDAFX_H_

#include <iostream>

#define SIMULATION_2D 1
//#define CONNECTED 1
//#define OBSTACLE 1
#define OUTPUT 1

#ifdef SIMULATION_2D
	#define _W 3 
	#define VISBLEW 600
	#define _H 2 //30
	#define GRIDSIZE 50 //15
	#define DIFFUSION 0.01
	#define VISCOSITY 0.01
	#define TIMESTEP 0.01
	#define ITERATION 10
	#define FRAMERATE 32
	#define DRAGSCALE 100

#ifdef CONNECTED
	#define IX(x, y) ( (x) == 0? _W + (y) * (_W+2) : ((x) == _W+1? 1 + (y) * (_W+2) : (x) + (y) * (_W+2)) )
#else	
	#define IX(x, y) ( (x) + (y) * (_W+2) )
#endif
	#define BOUNDED(x, y) ( (type[IX(int(x),int(y))] == SOLID || type[IX(int(x)+1,int(y)+1)] == SOLID)? false : true)
	
	//int IX(float x, float y);
#else
	#define _N 40
	#define GRIDSIZE 20
	#define DIFFUSION 0.0001
	#define VISCOSITY 0.0001
	#define TIMESTEP 0.01
	#define ITERATION 10
	#define FRAMERATE 32

	#define IX(x, y, z) ((x) + (y)*(_N+2) + (z)*(_N+2)*(_N+2) )
#endif

enum GRIDTYPE
{
	FLUID,
	AIR,
	SOLID
};

#define PI 3.14159265
#define LENGTH _N*GRIDSIZE
#define eps 1e-10
#define SWAP(x0, x) {float *tmp = x0; x0 = x; x = tmp;}

//system output defines
#define REPORT(X) std::cout << #X << ": " << X << std::endl
#define PRINT(X) std::cout << X << std::endl
#define PRINT_ONELINE(X) std::cout << X

#endif