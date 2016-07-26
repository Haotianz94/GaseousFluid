#include "stdafx.h"
#include "fluidCube.h"
#include <memory.h>
#include <math.h>
//#include "vl\VLfd.h"
#include "vl\VLf.h"
#include "vl\Solve.h"

#ifdef SIMULATION_2D

FluidCube2D::FluidCube2D(float diffusion, float viscosity, float dtime)
{
	size = (_W+2) * (_H+2);
	h = 1 / _H;
	h2 = h * h;
	dt = dtime;
	diff = diffusion;
	visc = viscosity;

	max_d = 0;
	max_vx = 0;
	max_vy = 0;

	d = new float [size];
	d0 = new float [size];
	Vx = new float [size]; 
	Vy = new float [size]; 
	Vx0 = new float [size]; 
	Vy0 = new float [size]; 
	type = new GRIDTYPE [size];

	//Advection using BFECC
	d_bf = new float [size];
	Vx_b = new float [size];
	Vx_f = new float [size];
	Vy_b = new float [size];
	Vy_f = new float [size];

	//Projection using Conjugate Gradient
	dir[0] = Pos(0, -1);
	dir[1] = Pos(-1, 0);
	dir[2] = Pos(1, 0);
	dir[3] = Pos(0, 1);
	pos2index = new int [size]; 
	neighNum = new int [size];
	neighbor = new int* [size];
	for(int i = 0; i < size; i ++)
		neighbor[i] = new int[4];

	for(int i = 0; i < size; i++)
	{
		pos2index[i] = -1;
		type[i] = SOLID;
	}
	//for(int y = 0; y <= _H+1; y++)
	//	type[IX(_W+1, y)] = AIR;	

#ifdef OBSTACLE
	int cx = 15;
	int cy = _H / 2;
	int R = _H * 0.25;
	fluidNum = 0;
	for(int y = 1; y <= _H; y++)
		for(int x = 1; x <= _W; x++)
		{
			float r = sqrtf((x-cx)*(x-cx) + (y-cy)*(y-cy));
			if(r < R)
				type[IX(x, y)] = SOLID;
			else
			{
				type[IX(x, y)] = FLUID;
				pos2index[IX(x, y)] = fluidNum ++;
			}
		}
	for(int y = 1; y <= _H; y++)
		for(int x = 1; x <= _W; x++)
			if(type[IX(x, y)] == SOLID)
			{
				for(int i = 0; i < 4; i++)
					if(type[IX(x+dir[i].x, y+dir[i].y)] == FLUID)
					{
						obstacle.push_back(Pos(x, y));
						break;
					}
			}
			else //FLUID
			{
				neighNum[IX(x, y)] = 0;
				for(int i = 0; i < 4; i++)
				{
					int xx = x+dir[i].x;
					int yy = y+dir[i].y;
					neighbor[IX(x,y)][i] = pos2index[IX(xx, yy)];
					if(type[IX(xx,yy)] != SOLID)
						neighNum[IX(x, y)] ++;
				}
			}
#else
	fluidNum = 0;
	for(int y = 1; y <= _H; y++)
		for(int x = 1; x <= _W; x++)
		{
			type[IX(x, y)] = FLUID;
			pos2index[IX(x, y)] = fluidNum ++;
		}
	for(int y = 1; y <= _H; y++)
		for(int x = 1; x <= _W; x++)
		{
			neighNum[IX(x, y)] = 0;
			for(int i = 0; i < 4; i++)
			{
				int xx = x+dir[i].x;
				int yy = y+dir[i].y;
				neighbor[IX(x,y)][i] = pos2index[IX(xx, yy)];
				if(type[IX(xx,yy)] != SOLID)
					neighNum[IX(x, y)] ++;
			}
		}
				
#endif

	//init Matrix
	A = new SparseMatf(fluidNum, fluidNum);
	int index = 0;
	for(int y = 1; y <= _H; y++)
		for(int x = 1; x <= _W; x++)
		{
			if(type[IX(x, y)] != FLUID)
				continue;

			//b[index] = h2 / dt * 0.5 * (Vx[IX(x+1,y)]-Vx[IX(x-1,y)] + Vy[IX(x,y+1)]-Vy[IX(x,y-1)]);

			SparseVecf A0(fluidNum);
			A0.Begin();
			for(int i = 0; i < 2; i++)
			{
				int neighid = neighbor[IX(x,y)][i];
				if(neighid != -1)
				{
					A0.AddNZElt(neighid, -1);
				}
			}
			A0.AddNZElt(index, neighNum[IX(x, y)]);
			for(int i = 2; i < 4; i++)
			{
				int neighid = neighbor[IX(x,y)][i];
				if(neighid != -1)
				{
					A0.AddNZElt(neighid, -1);
				}
			}
			A0.End();
			//cout<<A0<<endl;
			//A0.SetElts(index, neighs, VL_SV_END);
			(*A)[index++] = A0;
		}
	//cout<<*A<<endl;

	memset(d, 0, sizeof(float) * size);
	memset(Vx, 0, sizeof(float) * size);
	memset(Vy, 0, sizeof(float) * size);
	memset(d0, 0, sizeof(float) * size);
	memset(Vx0, 0, sizeof(float) * size);
	memset(Vy0, 0, sizeof(float) * size);

	memset(d_bf, 0, sizeof(float) * size);
	memset(Vx_b, 0, sizeof(float) * size);
	memset(Vx_f, 0, sizeof(float) * size);
	memset(Vy_b, 0, sizeof(float) * size);
	memset(Vy_f, 0, sizeof(float) * size);
}

FluidCube2D::~FluidCube2D()
{
	delete [] d;
	delete [] d0;
	delete [] Vx;
	delete [] Vy;
	delete [] Vx0;
	delete [] Vy0;
	delete [] type;

	delete [] pos2index;
	for(int i = 0; i < size; i++)
		delete[] neighbor[i];
	delete [] neighbor;
	delete [] neighNum;

	delete [] d_bf;
	delete [] Vx_b;
	delete [] Vx_f;
	delete [] Vy_b;
	delete [] Vy_f;

}

void FluidCube2D::dens_step(float *amount)
{
	addDensity(amount);
	
	SWAP(d0, d);
	diffuseDensity();
	
	//REPORT(Vx[IX(10,10)]);
	//REPORT(Vy[IX(10,10)]);

	SWAP(d0, d);
	advectDensity();
	
	//REPORT(d0[IX(10,10)]);
	//REPORT(d[IX(10,10)]);
}

void FluidCube2D::addDensity(float *amount)
{
	add_source(d, amount);
}

void FluidCube2D::diffuseDensity()
{
	diffuse(0, d0, d, diff);
}

void FluidCube2D::advectDensity()
{
	advect(0, d0, d, true);
	advect(0, d, d_bf, false);
	for(int i = 0; i < size; i++)
		d0[i] = d0[i] + (d0[i] - d_bf[i]) * 0.5;
	advect(0, d0, d, true);
}

void FluidCube2D::vel_step(float *amountX, float *amountY)
{
	addVelocity(amountX, amountY);
	
	SWAP(Vx0, Vx);
	SWAP(Vy0, Vy);
	diffuseVelosity();
	
	projectVelosity();

	SWAP(Vx0, Vx);
	SWAP(Vy0, Vy);
	advectVelosity();
	
	projectVelosity();
}

void FluidCube2D::addVelocity(float *amountX, float *amountY)
{
	add_source(Vx, amountX);
	add_source(Vy, amountY);
}

void FluidCube2D::diffuseVelosity()
{
	diffuse(1, Vx0, Vx, visc);
	diffuse(2, Vy0, Vy, visc);
}

void FluidCube2D::advectVelosity()
{
	advect(1, Vx0, Vx_b, true);
	advect(1, Vx_b, Vx_f, false);
	output(Vx0);
	output(Vx_f);
	for(int i = 0; i < size; i++)
		Vx0[i] = Vx0[i] + (Vx0[i] - Vx_f[i]) * 0.5;
	advect(1, Vx0, Vx, true);


	advect(2, Vy0, Vy_b, true);
	advect(2, Vy_b, Vy_f, false);
	for(int i = 0; i < size; i++)
		Vy0[i] = Vy0[i] + (Vy0[i] - Vy_f[i]) * 0.5;
	advect(2, Vy0, Vy, true);
}

void FluidCube2D::projectVelosity()
{

	max_vx = 0;
	max_vy = 0;
	float *div = Vy0;

#ifdef GAUSS_SEIDEL
	//Gauss_Seidel
	float *p = Vx0;
	//float *div = Vy0;

	for(int y = 1; y <= _H; y++)
		for(int x = 1; x <= _W; x++)
		{
			if(type[IX(x, y)] != FLUID)
				continue;

			div[IX(x, y)] = -0.5 * h * (Vx[IX(x+1,y)]-Vx[IX(x-1,y)] + Vy[IX(x,y+1)]-Vy[IX(x,y-1)]);
			p[IX(x, y)] = 0;
		}
	set_bnd(0, div);
	set_bnd(0, p);

	//output(div);
	
	for(int k = 0; k < ITERATION; k++)
	{
		for(int y = 1; y <= _H; y++)
			for(int x = 1; x <= _W; x++)
				if(type[IX(x, y)] == FLUID)
					p[IX(x, y)] = (div[IX(x,y)] + p[IX(x-1,y)] + p[IX(x+1,y)] + p[IX(x,y-1)] + p[IX(x,y+1)]) / 4;
		//keep balance
		/*if(k % 1 == 0)
		{
   			for(int y = 1; y <= _H; y++)
				for(int x = 1; x <= _W; x++)
					p[IX(x, y)] = (div[IX(x,y)] + p[IX(x-1,y)] + p[IX(x+1,y)] + p[IX(x,y-1)] + p[IX(x,y+1)]) / 4;
		}
		else if(k % 4 == 1)
		{
			for(int y = _N; y >= 1; y--)
				for(int x = _N; x >= 1; x--)
					p[IX(x, y)] = (div[IX(x,y)] + p[IX(x-1,y)] + p[IX(x+1,y)] + p[IX(x,y-1)] + p[IX(x,y+1)]) / 4;
		}
		else if(k % 4 == 2)
		{
			for(int y = 1; y <= _N; y++)
				for(int x = _N; x >= 1; x--)
					p[IX(x, y)] = (div[IX(x,y)] + p[IX(x-1,y)] + p[IX(x+1,y)] + p[IX(x,y-1)] + p[IX(x,y+1)]) / 4;
		}
		else
		{
			for(int y = _N; y >= 1; y--)
				for(int x = 1; x <= _N; x++)
					p[IX(x, y)] = (div[IX(x,y)] + p[IX(x-1,y)] + p[IX(x+1,y)] + p[IX(x,y-1)] + p[IX(x,y+1)]) / 4;
		}*/
		set_bnd(0, p);
	}
 
	for(int y = 1; y <= _H; y++)
		for(int x = 1; x <= _W; x++)
		{
			if(type[IX(x, y)] != FLUID)
				continue;

			Vx[IX(x, y)] -= 0.5 * (p[IX(x+1,y)] - p[IX(x-1,y)]) / h;
			Vy[IX(x, y)] -= 0.5 * (p[IX(x,y+1)] - p[IX(x,y-1)]) / h;

			if(fabsf(Vx[IX(x, y)]) > max_vx)
				max_vx = fabsf(Vx[IX(x, y)]);
			if(fabsf(Vy[IX(x, y)]) > max_vy)
				max_vy = fabsf(Vy[IX(x, y)]);
		}
	set_bnd(1, Vx);
	set_bnd(2, Vy);
#else
	//Conjugate Gradient

	//check div before project
	/*for(int y = 1; y <= _H; y++)
		for(int x = 1; x <= _W; x++)
		{
			if(type[IX(x, y)] != FLUID)
				continue;
			div[IX(x, y)] = 0.5 * (Vx[IX(x+1,y)]-Vx[IX(x-1,y)] + Vy[IX(x,y+1)]-Vy[IX(x,y-1)]);
		}
	output(div);*/

	Vecf b(fluidNum);
	int index = 0;
	for(int y = 1; y <= _H; y++)
		for(int x = 1; x <= _W; x++)
		{
			if(type[IX(x, y)] != FLUID)
				continue;

			b[index] =  -0.5 * h * (Vx[IX(x+1,y)]-Vx[IX(x-1,y)] + Vy[IX(x,y+1)]-Vy[IX(x,y-1)]);

			SparseVecf A0(fluidNum);
			A0.Begin();
			for(int i = 0; i < 2; i++)
			{
				int neighid = neighbor[IX(x,y)][i];
				if(neighid != -1)
				{
					A0.AddNZElt(neighid, -1);
				}
			}
			A0.AddNZElt(index, neighNum[IX(x, y)]);
			for(int i = 2; i < 4; i++)
			{
				int neighid = neighbor[IX(x,y)][i];
				if(neighid != -1)
				{
					A0.AddNZElt(neighid, -1);
				}
			}
			A0.End();
			//cout<<A0<<endl;
			//A0.SetElts(index, neighs, VL_SV_END);
			(*A)[index++] = A0;
		}
	Vecf p(b);
	SolveConjGrad(*A, p, b, 1e-7);
	
	//for(int i = 100; i < 110; i++)
	//	cout<<p[i]<<' ';
	//cout<<endl;
	//cout<<p<<endl;

	for(int y = 1; y <= _H; y++)
		for(int x = 1; x <= _W; x++)
		{
			if(type[IX(x, y)] != FLUID)
				continue;

			float p1, p2;
			if(type[IX(x+1, y)] == AIR)
				p2 = 0;
			else if(type[IX(x+1, y)] == FLUID)
				p2 = p[pos2index[IX(x+1,y)]];
			else
				p2 = p[pos2index[IX(x,y)]];
			if(type[IX(x-1, y)] == AIR)
				p1 = 0;
			else if(type[IX(x-1, y)] == FLUID)
				p1 = p[pos2index[IX(x-1,y)]];
			else
				p1 = p[pos2index[IX(x,y)]];
			Vx[IX(x, y)] -= 0.5 * (p2 - p1) / h;
			
			if(type[IX(x, y+1)] == AIR)
				p2 = 0;
			else if(type[IX(x, y+1)] == FLUID)
				p2 = p[pos2index[IX(x,y+1)]];
			else
				p2 = p[pos2index[IX(x,y)]];
			if(type[IX(x, y-1)] == AIR)
				p1 = 0;
			else if(type[IX(x, y-1)] == FLUID)
				p1 = p[pos2index[IX(x,y-1)]];
			else
				p1 = p[pos2index[IX(x,y)]];
			Vy[IX(x, y)] -= 0.5 * (p2 - p1) / h;

			if(fabsf(Vx[IX(x, y)]) > max_vx)
				max_vx = fabsf(Vx[IX(x, y)]);
			if(fabsf(Vy[IX(x, y)]) > max_vy)
				max_vy = fabsf(Vy[IX(x, y)]);
		}
	set_bnd(1, Vx);
	set_bnd(2, Vy);


	//check div after project
	/*for(int y = 1; y <= _H; y++)
		for(int x = 1; x <= _W; x++)
		{
			if(type[IX(x, y)] != FLUID)
				continue;
			div[IX(x, y)] = 0.5 * (Vx[IX(x+1,y)]-Vx[IX(x-1,y)] + Vy[IX(x,y+1)]-Vy[IX(x,y-1)]);
		}
	output(div);
	int stop = 0;*/
#endif

}

void FluidCube2D::add_source(float *x, float *amount)
{
	for(int i = 0; i < size; i++)
		x[i] += dt*amount[i];
}

void FluidCube2D::diffuse(int b, float *u0, float *u, float diffusion)
{
	float a = dt * diffusion / h2;
	for(int k = 0; k < ITERATION; k++)
	{
		for(int y = 1; y <= _H; y++)
			for(int x = 1; x <= _W; x++)
				if(type[IX(x, y)] == FLUID)
					u[IX(x, y)] = (u0[IX(x, y)] + a * (u[IX(x-1, y)]+u[IX(x+1, y)]+u[IX(x,y-1)]+u[IX(x,y+1)])) / (1+4*a);
		//keep balance
		/*if(k % 4 == 0)
		{
   			for(int y = 1; y <= _H; y++)
				for(int x = 1; x <= _W; x++)
				{
					u[IX(x, y)] = (u0[IX(x, y)] + a * (u[IX(x-1, y)]+u[IX(x+1, y)]+u[IX(x,y-1)]+u[IX(x,y+1)])) / (1+4*a);
				}
		}
		else if(k % 4 == 1)
		{
			for(int y = _N; y >= 1; y--)
				for(int x = _N; x >= 1; x--)
				{
					u[IX(x, y)] = (u0[IX(x, y)] + a * (u[IX(x-1, y)]+u[IX(x+1, y)]+u[IX(x,y-1)]+u[IX(x,y+1)])) / (1+4*a);
				}
		}
		else if(k % 4 == 2)
		{
			for(int y = 1; y <= _N; y++)
				for(int x = _N; x >= 1; x--)
				{
					u[IX(x, y)] = (u0[IX(x, y)] + a * (u[IX(x-1, y)]+u[IX(x+1, y)]+u[IX(x,y-1)]+u[IX(x,y+1)])) / (1+4*a);
				}
		}
		else
		{
			for(int y = _N; y >= 1; y--)
				for(int x = 1; x <= _N; x++)
				{
					u[IX(x, y)] = (u0[IX(x, y)] + a * (u[IX(x-1, y)]+u[IX(x+1, y)]+u[IX(x,y-1)]+u[IX(x,y+1)])) / (1+4*a);
				}
		}*/
		set_bnd(b, u);
	}
}

void FluidCube2D::advect(int b, float *u0, float *u, bool backward)
{
	max_d = 0;
	float dt0 = dt / h;
	for(int y = 1; y <= _H; y++)
		for(int x = 1; x <= _W; x++)
		{
			if(type[IX(x, y)] != FLUID)
				continue;

			float x0, y0;
			if(backward)
			{
				x0 = x - dt0*Vx[IX(x,y)];
				y0 = y - dt0*Vy[IX(x,y)];
			}
			else
			{
				x0 = x + dt0*Vx[IX(x,y)];
				y0 = y + dt0*Vy[IX(x,y)];
			}

#ifdef CONNECTED
			if(x0 <= 0)
				x0 += _W;
			else if(x0 >= _W+1)
				x0 -= _W;
#else
			if(x0 < 0.5)
				x0 = 0.5;
			else if(x0 > _W + 0.5)
				x0 = _W + 0.5;
#endif
			if(y0 < 0.5)
				y0 = 0.5;
			else if(y0 > _H + 0.5)
				y0 = _H + 0.5;

			//if(!BOUNDED(x0, y0))
			//	continue;

			int i0 = int(x0), i1 = i0 + 1;
			int j0 = int(y0), j1 = j0 + 1;
			float s1 = x0 - i0, s0 = 1 - s1;
			float t1 = y0 - j0, t0 = 1 - t1;

			u[IX(x,y)] = s0 * (t0*u0[IX(i0,j0)] + t1*u0[IX(i0,j1)]) +
						 s1 * (t0*u0[IX(i1,j0)] + t1*u0[IX(i1,j1)]);

			if(b == 0)
			{
				if(u[IX(x,y)] == 0)
					continue;

				//float v = log(u[IX(x,y)]);
				float v = u[IX(x,y)];
				if(v > max_d)
					max_d = v;
			}
		}

	set_bnd(b, u);
}

void FluidCube2D::set_bnd(int b, float *x)
{
	for(unsigned i = 0; i < obstacle.size(); i++)
	{
		int ox = obstacle[i].x;
		int oy = obstacle[i].y;
		if(b == 0)
		{
			float td = 0;
			int nei = 0;
			for(int j = 0; j < 4; j++)
			{
				int xx = ox + dir[j].x;
				int yy = oy + dir[j].y;
				if(type[IX(xx,yy)] == FLUID)
				{
					td += x[IX(xx,yy)];
					nei++;
				}
			}
			x[IX(ox, oy)] = td / nei;
		}
		else if(b == 1)
		{
			if(type[IX(ox-1, oy)] == FLUID)
				x[IX(ox, oy)] = -x[IX(ox-1, oy)];
			else if(type[IX(ox+1, oy)] == FLUID)
				x[IX(ox, oy)] = -x[IX(ox+1, oy)];
			else if(type[IX(ox, oy-1)] == FLUID)
				x[IX(ox, oy)] = x[IX(ox-1, oy)];
			else
				x[IX(ox, oy)] = x[IX(ox+1, oy)];
		}
		else
		{
			if(type[IX(ox, oy-1)] == FLUID)
				x[IX(ox, oy)] = -x[IX(ox, oy-1)];
			else if(type[IX(ox, oy+1)] == FLUID)
				x[IX(ox, oy)] = -x[IX(ox, oy+1)];
			else if(type[IX(ox-1, oy)] == FLUID)
				x[IX(ox, oy)] = x[IX(ox-1, oy)];
			else
				x[IX(ox, oy)] = x[IX(ox+1, oy)];
		}
	}

	for(int i = 1; i <= _W; i++)
	{
		x[IX(i, 0)] = b==2? -x[IX(i,1)] : x[IX(i,1)];
		x[IX(i, _H+1)] = b==2? -x[IX(i,_H)] : x[IX(i,_H)];
 	}

#ifndef CONNECTED
	
	for(int i = 1; i <= _H; i++)
	{
		x[IX(0, i)] = b==1? -x[IX(1,i)] : x[IX(1,i)];
		//x[IX(_W+1, i)] = b==1? -x[IX(_W,i)] : x[IX(_W,i)];
 	}
	x[IX(0, 0)] = 0.5 * (x[IX(1, 0)] + x[IX(0, 1)]);
	x[IX(0, _H+1)] = 0.5 * (x[IX(1, _H+1)] + x[IX(0, _H)]);
	x[IX(_W+1, 0)] = 0.5 * (x[IX(_W, 0)] + x[IX(_W+1, 1)]);
	x[IX(_W+1, _H+1)] = 0.5 * (x[IX(_W, _H+1)] + x[IX(_W+1, _H)]);

#endif
}

void FluidCube2D::simulate(bool idle)
{
	if(idle)
	{
		memset(d0, 0, sizeof(float) * size);
		memset(Vx0, 0, sizeof(float) * size);
		memset(Vy0, 0, sizeof(float) * size);
	}

	vel_step(Vx0, Vy0);
	dens_step(d0);

	draw_dens();
}

void FluidCube2D::output(float *u)
{
#ifdef OUTPUT 
	for(int y = 10; y <= 15; y++)
		for(int x = 40; x <= 45; x++)
		{
			std::cout<<u[IX(x, y)]<<' ';
			if(x == 15)
				std::cout<<std::endl;
		}
	PRINT("=================================\n");
#endif
}

void FluidCube2D::draw_dens()
{
	//identify that we are currently modifying the projection matrix
	glMatrixMode(GL_PROJECTION);
	//load the identity matrix (clear old matrix)
	glLoadIdentity();

	//get the current window width/height
	int w = glutGet(GLUT_WINDOW_WIDTH);
	int h = glutGet(GLUT_WINDOW_HEIGHT);

	//initialize the screen coordinates
	//glOrtho(-w/2, w/2-1, -h/2, h/2-1, -1, 1);
	glOrtho(-10, w+10, -10, h+10, -1, 1);

	//now we are editing the modelview matrix
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	//clear the screen to a desired color in range [0-1] RGBA
	glClearColor(0.5, 0.5, 0.5, 1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	//enable blending for translucency
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_LINE_SMOOTH);

	REPORT(max_d);
	REPORT(max_vx);
	REPORT(max_vy);

	//
	//max_d *= 2.0/3;

	int W = VISBLEW/GRIDSIZE;
	for(int i = 0; i < W; i++)
		for(int j = 0; j < _H; j++)
		{
			int x = i + 1;
			int y = j + 1;
			float color;

			if(type[IX(x, y)] == SOLID)
				glColor3f(0, 0.5, 0); 
			else if(d[IX(x, y)] == 0)
				glColor3f(0, 0, 0);
			else
			{
				//color = (log(d[IX(x, (h-y))]) - min) / gap;
				color = d[IX(x, y)] / max_d;
				glColor3f(color, color, color);
			}
			glBegin(GL_QUADS);
			glVertex2f(i*GRIDSIZE, j*GRIDSIZE);
			glVertex2f((i+1)*GRIDSIZE, j*GRIDSIZE);
			glVertex2f((i+1)*GRIDSIZE, (j+1)*GRIDSIZE);
			glVertex2f(i*GRIDSIZE, (j+1)*GRIDSIZE);
			glEnd();

			if(GRIDSIZE >= 10 && type[IX(x, y)] == FLUID)
				draw_velo(i, j, Vx[IX(x, y)], Vy[IX(x, y)]);
		}

	glutSwapBuffers();
}

void FluidCube2D::draw_velo(int i, int j, float vx, float vy)
{
	float cx = (i+0.5) * GRIDSIZE;
	float cy = (j+0.5) * GRIDSIZE;
	float vl = sqrtf(vx*vx + vy*vy);
	//float seita = atanf(fabs(vy / vx));

	float max_v = max_vx > max_vy? max_vx : max_vy;
	if(max_v == 0)
		return;

	/*if(fabs(vx) < eps && fabs(vy) < eps)
		return;
	int LOG = 10;
	max_v = logf(max_v) + LOG;
	float dx = 0.6 * GRIDSIZE * (logf(fabsf(vx))+LOG) / max_v * ((vx>0)? 1:-1);
	float dy = 0.6 * GRIDSIZE * (logf(fabsf(vy))+LOG) / max_v * ((vy>0)? 1:-1);*/
	
	float dx = 0.5 * GRIDSIZE * vx / vl;
	float dy = 0.5 * GRIDSIZE * vy / vl;

	//float dx = 20 * GRIDSIZE * vx ;
	//float dy = 20 * GRIDSIZE * vy ;

	glBegin(GL_LINES);
	glColor3f(1.0f, 0.0f, 0.0f);
	glVertex2f(cx, cy);
	//glColor3f(0.0f, 0.0f, 1.0f);
	glVertex2f(cx + dx, cy + dy);
	
	glColor3f(0.0f, 0.0f, 1.0f);
	glVertex2f(cx, cy);
	glVertex2f(cx - vx/vl*GRIDSIZE*0.2, cy - vy/vl*GRIDSIZE*0.2);
	glEnd();
}

#endif