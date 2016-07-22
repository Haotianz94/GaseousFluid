#include "stdafx.h"
#include "fluidCube.h"
#include <memory.h>
#include <math.h>

#ifndef SIMULATION_2D

extern float px;
extern float py;
extern float pz;

FluidCube3D::FluidCube3D(float diffusion, float viscosity, float dtime)
{
	size = (_N+2) * (_N+2) * (_N+2);
	h = 1.0 / _N;
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
	Vz = new float [size]; 
	Vx0 = new float [size]; 
	Vy0 = new float [size]; 
	Vz0 = new float [size];

	memset(d, 0, sizeof(float) * size);
	memset(Vx, 0, sizeof(float) * size);
	memset(Vy, 0, sizeof(float) * size);
	memset(Vz, 0, sizeof(float) * size);
	memset(d0, 0, sizeof(float) * size);
	memset(Vx0, 0, sizeof(float) * size);
	memset(Vy0, 0, sizeof(float) * size);
	memset(Vz0, 0, sizeof(float) * size);
}

FluidCube3D::~FluidCube3D()
{
	delete [] d;
	delete [] d0;
	delete [] Vx;
	delete [] Vy;
	delete [] Vz;
	delete [] Vx0;
	delete [] Vy0;
	delete [] Vz0;
}

void FluidCube3D::dens_step(float *amount)
{
	addDensity(amount);
	
	SWAP(d0, d);
	diffuseDensity();
	
	REPORT(Vx[IX(5,5,5)]);
	REPORT(Vy[IX(5,5,5)]);
	REPORT(Vz[IX(5,5,5)]);

	SWAP(d0, d);
	advectDensity();
	
	REPORT(d0[IX(5,5,5)]);
	REPORT(d[IX(5,5,5)]);
}

void FluidCube3D::addDensity(float *amount)
{
	add_source(d, amount);
}

void FluidCube3D::diffuseDensity()
{
	diffuse(0, d0, d, diff);
}

void FluidCube3D::advectDensity()
{
	advect(0, d0, d);
}

void FluidCube3D::vel_step(float *amountX, float *amountY, float *amountZ)
{
	addVelocity(amountX, amountY, amountZ);

	SWAP(Vx0, Vx);
	SWAP(Vy0, Vy);
	SWAP(Vz0, Vz);
	diffuseVelosity();
	
	projectVelosity();

	SWAP(Vx0, Vx);
	SWAP(Vy0, Vy);
	SWAP(Vz0, Vz);
	advectVelosity();

	projectVelosity();
}

void FluidCube3D::addVelocity(float *amountX, float *amountY, float *amountZ)
{
	add_source(Vx, amountX);
	add_source(Vy, amountY);
	add_source(Vz, amountZ);
}

void FluidCube3D::diffuseVelosity()
{
	diffuse(1, Vx0, Vx, visc);
	diffuse(2, Vy0, Vy, visc);
	diffuse(3, Vz0, Vz, visc);
}

void FluidCube3D::advectVelosity()
{
	advect(1, Vx0, Vx);
	advect(2, Vy0, Vy);
	advect(3, Vz0, Vz);
}

void FluidCube3D::projectVelosity()
{
	max_vx = 0;
	max_vy = 0;
	max_vz = 0;

	float *p = Vx0;
	float *div = Vy0;

	for(int z = 1; z <= _N; z++)
		for(int y = 1; y <= _N; y++)
			for(int x = 1; x <= _N; x++)
			{
				div[IX(x, y, z)] = -0.5 * h * (Vx[IX(x+1,y,z)]-Vx[IX(x-1,y,z)] + Vy[IX(x,y+1,z)]-Vy[IX(x,y-1,z)] +
					Vz[IX(x,y,z+1)] - Vz[IX(x,y,z-1)]);
				p[IX(x, y, z)] = 0;
			}
	set_bnd(0, div);
	set_bnd(0, p);
	
	for(int k = 0; k < ITERATION; k++)
	{
		//if(k % 8 == 0)
		//{
		for(int z = 1; z <= _N; z++)
   			for(int y = 1; y <= _N; y++)
				for(int x = 1; x <= _N; x++)
					p[IX(x, y, z)] = (div[IX(x,y,z)] + p[IX(x-1,y,z)] + p[IX(x+1,y,z)] + p[IX(x,y-1,z)] + p[IX(x,y+1,z)]
						+ p[IX(x,y,z-1)] + p[IX(x,y,z+1)]) / 6;
		//}
		set_bnd(0, p);
	}

	for(int z = 1; z <= _N; z++)
		for(int y = 1; y <= _N; y++)
			for(int x = 1; x <= _N; x++)
			{
				Vx[IX(x, y, z)] -= 0.5 * (p[IX(x+1,y,z)] - p[IX(x-1,y,z)]) / h;
				Vy[IX(x, y, z)] -= 0.5 * (p[IX(x,y+1,z)] - p[IX(x,y-1,z)]) / h;
				Vz[IX(x, y, z)] -= 0.5 * (p[IX(x,y,z+1)] - p[IX(x,y,z-1)]) / h;

				if(fabsf(Vx[IX(x, y, z)]) > max_vx)
					max_vx = fabsf(Vx[IX(x, y, z)]);
				if(fabsf(Vy[IX(x, y, z)]) > max_vy)
					max_vy = fabsf(Vy[IX(x, y, z)]);
				if(fabsf(Vz[IX(x, y, z)]) > max_vz)
					max_vz = fabsf(Vz[IX(x, y, z)]);
			}
	set_bnd(1, Vx);
	set_bnd(2, Vy);
	set_bnd(3, Vz);
}

void FluidCube3D::add_source(float *x, float *amount)
{
	for(int i = 0; i < size; i++)
		x[i] += dt*amount[i];
}

void FluidCube3D::diffuse(int b, float *u0, float *u, float diffusion)
{
	float a = dt * diffusion * _N * _N * _N;
	for(int k = 0; k < ITERATION; k++)
	{
		//if(k % 8 == 0)
		//{
		for(int z = 1; z <= _N; z++)
   			for(int y = 1; y <= _N; y++)
				for(int x = 1; x <= _N; x++)
					u[IX(x, y, z)] = (u0[IX(x, y, z)] + a * (u[IX(x-1, y, z)]+u[IX(x+1, y, z)]+u[IX(x,y-1,z)]+
						u[IX(x,y+1,z)]+u[IX(x,y,z-1)]+u[IX(x,y,z+1)])) / (1+6*a);
		//}
		set_bnd(b, u);
	}
}

void FluidCube3D::advect(int b, float *u0, float *u)
{
	max_d = 0;
	float dt0 = dt / h;
	for(int z = 1; z <= _N; z++)
		for(int y = 1; y <= _N; y++)
			for(int x = 1; x <= _N; x++)
			{
				float x0 = x - dt0*Vx[IX(x,y,z)];
				float y0 = y - dt0*Vy[IX(x,y,z)];
				float z0 = z - dt0*Vz[IX(x,y,z)];

				if(x0 < 0.5)
					x0 = 0.5;
				else if(x0 > _N + 0.5)
					x0 = _N + 0.5;
				if(y0 < 0.5)
					y0 = 0.5;
				else if(y0 > _N + 0.5)
					y0 = _N + 0.5;
				if(z0 < 0.5)
					z0 = 0.5;
				else if(z0 > _N + 0.5)
					z0 = _N + 0.5;

				int i0 = int(x0), i1 = i0 + 1;
				int j0 = int(y0), j1 = j0 + 1;
				int k0 = int(z0), k1 = k0 + 1;
				float s1 = x0 - i0, s0 = 1 - s1;
				float t1 = y0 - j0, t0 = 1 - t1;
				float r1 = z0 - k0, r0 = 1 - r1;

				float a1 = r0 * u0[IX(i0, j0, k0)] + r1 * u0[IX(i0, j0, k1)];
				float a2 = r0 * u0[IX(i0, j1, k0)] + r1 * u0[IX(i0, j1, k1)];
				float b1 = r0 * u0[IX(i1, j0, k0)] + r1 * u0[IX(i1, j0, k1)];
				float b2 = r0 * u0[IX(i1, j1, k0)] + r1 * u0[IX(i1, j1, k1)];

				float c1 = t0 * a1 + t1 * a2;
				float c2 = t0 * b1 + t1 * b2;

				u[IX(x, y, z)] = s0 * c1 + s1 * c2;

				if(b == 0)
				{
					if(u[IX(x,y,z)] == 0)
						continue;

					//float v = log(u[IX(x,y)]);
					float v = u[IX(x,y,z)];
					if(v > max_d)
						max_d = v;
				}
			}
	set_bnd(b, u);
}

void FluidCube3D::set_bnd(int b, float *x)
{
	// 6 face
	for(int i = 1; i <= _N; i++)
		for(int j = 1; j <= _N; j++)
		{
			x[IX(0, i, j)] = b==1? -x[IX(1,i,j)] : x[IX(1,i,j)];
			x[IX(_N+1, i, j)] = b==1? -x[IX(_N,i,j)] : x[IX(_N,i,j)];

			x[IX(i, 0, j)] = b==2? -x[IX(i,1,j)] : x[IX(i,1,j)];
			x[IX(i, _N+1, j)] = b==2? -x[IX(i,_N,j)] : x[IX(i,_N,j)];

			x[IX(i, j, 0)] = b==3? -x[IX(i,j,1)] : x[IX(i,j,1)];
			x[IX(i, j, _N+1)] = b==3? -x[IX(i,j,_N)] : x[IX(i,j,_N)];
 		}
	// 12 edges
	for(int i = 1; i <= _N; i++)
	{
		x[IX(0, 0, i)] = 0.5 * (x[IX(0, 1, i)] + x[IX(1, 0, i)]);
		x[IX(0, _N+1, i)] = 0.5 * (x[IX(1, _N+1, i)] + x[IX(0, _N, i)]);
		x[IX(_N+1, 0, i)] = 0.5 * (x[IX(_N+1, 1, i)] + x[IX(_N, 0, i)]);
		x[IX(_N+1, _N+1, i)] = 0.5 * (x[IX(_N, _N+1, i)] + x[IX(_N+1, _N, i)]);

		x[IX(0, i, 0)] = 0.5 * (x[IX(0, i, 1)] + x[IX(1, i, 0)]);
		x[IX(0, i, _N+1)] = 0.5 * (x[IX(1, i, _N+1)] + x[IX(0, i, _N)]);
		x[IX(_N+1, i, 0)] = 0.5 * (x[IX(_N+1, i, 1)] + x[IX(_N, i, 0)]);
		x[IX(_N+1, i, _N+1)] = 0.5 * (x[IX(_N, i, _N+1)] + x[IX(_N+1, i, _N)]);

		x[IX(i, 0, 0)] = 0.5 * (x[IX(i, 0, 1)] + x[IX(i, 1, 0)]);
		x[IX(i, 0, _N+1)] = 0.5 * (x[IX(i, 1, _N+1)] + x[IX(i, 0, _N)]);
		x[IX(i, _N+1, 0)] = 0.5 * (x[IX(i, _N+1, 1)] + x[IX(i, _N, 0)]);
		x[IX(i, _N+1, _N+1)] = 0.5 * (x[IX(i, _N, _N+1)] + x[IX(i, _N+1, _N)]);
	}
	// 8 vertexs
	x[IX(0, 0, 0)]		 = 0.33f * (x[IX(1, 0, 0)]
                                  + x[IX(0, 1, 0)]
                                  + x[IX(0, 0, 1)]);
    x[IX(0, _N+1, 0)]     = 0.33f * (x[IX(1, _N+1, 0)]
                                  + x[IX(0, _N, 0)]
                                  + x[IX(0, _N+1, 1)]);
    x[IX(0, 0, _N+1)]     = 0.33f * (x[IX(1, 0, _N+1)]
                                  + x[IX(0, 1, _N+1)]
                                  + x[IX(0, 0, _N)]);
    x[IX(0, _N+1, _N+1)]   = 0.33f * (x[IX(1, _N+1, _N+1)]
                                  + x[IX(0, _N, _N+1)]
                                  + x[IX(0, _N+1, _N)]);
    x[IX(_N+1, 0, 0)]     = 0.33f * (x[IX(_N, 0, 0)]
                                  + x[IX(_N+1, 1, 0)]
                                  + x[IX(_N+1, 0, 1)]);
    x[IX(_N+1, _N+1, 0)]   = 0.33f * (x[IX(_N, _N+1, 0)]
                                  + x[IX(_N+1, _N, 0)]
                                  + x[IX(_N+1, _N+1, 1)]);
    x[IX(_N+1, 0, _N+1)]   = 0.33f * (x[IX(_N, 0, _N+1)]
                                  + x[IX(_N+1, 1, _N+1)]
                                  + x[IX(_N+1, 0, _N)]);
    x[IX(_N+1, _N+1, _N+1)] = 0.33f * (x[IX(_N, _N+1, _N+1)]
                                  + x[IX(_N+1, _N, _N+1)]
                                  + x[IX(_N+1, _N+1, _N)]);
}

void FluidCube3D::simulate(bool idle)
{
	if(idle)
	{
		memset(d0, 0, sizeof(float) * size);
		memset(Vx0, 0, sizeof(float) * size);
		memset(Vy0, 0, sizeof(float) * size);
		memset(Vz0, 0, sizeof(float) * size);
	}

	vel_step(Vx0, Vy0, Vz0);
	dens_step(d0);

	draw_dens();
}

void FluidCube3D::output(float *u)
{
#ifdef OUTPUT 
	for(int y = 0; y <= _N+1; y++)
	{
		for(int x = 0; x <= _N+1; x++)
			printf("%.3f ", u[IX(x, y)]);
		printf("\n");
	}
	printf("=================================\n");
#endif
}

void FluidCube3D::draw_dens()
{
	glMatrixMode(GL_MODELVIEW);
	// Reset transformations
	glLoadIdentity();
	glScalef(0.02f, 0.02f, 0.02f);

	// Set the camera
	gluLookAt(	px, py, pz,
				0, 0, 0,
				0.0f, 1.0f, 0.0f);

	//clear the screen to a desired color in range [0-1] RGBA
	glClearColor(0.5, 0.5, 0.5, 1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	//enable blending for translucency
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_LINE_SMOOTH);

	REPORT(max_d);

	glTranslatef(-LENGTH/2, -LENGTH/2, -LENGTH/2);
	for(int k = 0; k < _N; k++)
		for(int j = 0; j < _N; j++)
			for(int i = 0; i < _N; i++)
			{
				int x = i + 1;
				int y = j + 1;
				int z = k + 1;

				float color;
				if(fabsf(d0[IX(x, y, z)]) == 0)
					continue;
				else
					//color = (log(d[IX(x, (h-y))]) - min) / gap;
					color = d0[IX(x, y, z)] / max_d;
				//PRINT(d0[IX(i,j,k)]);
				glColor4f(0, 0, 0.5, color);
				
				//draw cube 
				glBegin(GL_QUADS);
				//hold k
				glVertex3f(i*GRIDSIZE, j*GRIDSIZE, k*GRIDSIZE);
				glVertex3f(i*GRIDSIZE, (j+1)*GRIDSIZE, k*GRIDSIZE);
				glVertex3f((i+1)*GRIDSIZE, (j+1)*GRIDSIZE, k*GRIDSIZE);
				glVertex3f((i+1)*GRIDSIZE, j*GRIDSIZE, k*GRIDSIZE);
				//hold k+1
				glVertex3f(i*GRIDSIZE, j*GRIDSIZE, (k+1)*GRIDSIZE);
				glVertex3f((i+1)*GRIDSIZE, j*GRIDSIZE, (k+1)*GRIDSIZE);
				glVertex3f((i+1)*GRIDSIZE, (j+1)*GRIDSIZE, (k+1)*GRIDSIZE);
				glVertex3f(i*GRIDSIZE, (j+1)*GRIDSIZE, (k+1)*GRIDSIZE);
				//hold j
				glVertex3f(i*GRIDSIZE, j*GRIDSIZE, k*GRIDSIZE);
				glVertex3f((i+1)*GRIDSIZE, j*GRIDSIZE, k*GRIDSIZE);
				glVertex3f((i+1)*GRIDSIZE, j*GRIDSIZE, (k+1)*GRIDSIZE);
				glVertex3f(i*GRIDSIZE, j*GRIDSIZE, (k+1)*GRIDSIZE);
				//hold j+1
				glVertex3f(i*GRIDSIZE, (j+1)*GRIDSIZE, k*GRIDSIZE);
				glVertex3f(i*GRIDSIZE, (j+1)*GRIDSIZE, (k+1)*GRIDSIZE);
				glVertex3f((i+1)*GRIDSIZE, (j+1)*GRIDSIZE, (k+1)*GRIDSIZE);
				glVertex3f((i+1)*GRIDSIZE, (j+1)*GRIDSIZE, k*GRIDSIZE);
				//hold i
				glVertex3f(i*GRIDSIZE, j*GRIDSIZE, k*GRIDSIZE);
				glVertex3f(i*GRIDSIZE, j*GRIDSIZE, (k+1)*GRIDSIZE);
				glVertex3f(i*GRIDSIZE, (j+1)*GRIDSIZE, (k+1)*GRIDSIZE);
				glVertex3f(i*GRIDSIZE, (j+1)*GRIDSIZE, k*GRIDSIZE);
				//hold i+1
				glVertex3f((i+1)*GRIDSIZE, j*GRIDSIZE, k*GRIDSIZE);
				glVertex3f((i+1)*GRIDSIZE, (j+1)*GRIDSIZE, k*GRIDSIZE);
				glVertex3f((i+1)*GRIDSIZE, (j+1)*GRIDSIZE, (k+1)*GRIDSIZE);
				glVertex3f((i+1)*GRIDSIZE, j*GRIDSIZE, (k+1)*GRIDSIZE);
				glEnd();

				//if(GRIDSIZE >= 10)
				//	draw_velo(i, j, k, Vx[IX(x, y, z)], Vy[IX(x, y, z)], Vz[IX(x, y, z)]);

			}
	//draw_grid();

	//draw box
	glColor3f(1, 1, 1);
	glBegin(GL_LINES);
	glVertex3f(0, 0, 0);
	glVertex3f(0, LENGTH, 0);
	glVertex3f(0, LENGTH, 0);
	glVertex3f(LENGTH, LENGTH, 0);
	glVertex3f(LENGTH, LENGTH, 0);
	glVertex3f(LENGTH, 0, 0);
	glVertex3f(LENGTH, 0, 0);
	glVertex3f(0, 0, 0);

	glVertex3f(0, 0, LENGTH);
	glVertex3f(0, LENGTH, LENGTH);
	glVertex3f(0, LENGTH, LENGTH);
	glVertex3f(LENGTH, LENGTH, LENGTH);
	glVertex3f(LENGTH, LENGTH, LENGTH);
	glVertex3f(LENGTH, 0, LENGTH);
	glVertex3f(LENGTH, 0, LENGTH);
	glVertex3f(0, 0, LENGTH);

	glVertex3f(0, 0, 0);
	glVertex3f(0, 0, LENGTH);
	glVertex3f(0, LENGTH, LENGTH);
	glVertex3f(0, LENGTH, 0);
	glVertex3f(LENGTH, LENGTH, LENGTH);
	glVertex3f(LENGTH, LENGTH, 0);
	glVertex3f(LENGTH, 0, LENGTH);
	glVertex3f(LENGTH, 0, 0);
	glEnd();

	glutSwapBuffers();
}

void FluidCube3D::draw_velo(int i, int j, int k, float vx, float vy, float vz)
{
	float cx = (i+0.5) * GRIDSIZE;
	float cy = (j+0.5) * GRIDSIZE;
	float cz = (k+0.5) * GRIDSIZE;

	/*float vl = sqrtf(vx*vx + vy*vy +vz*vz);
	float fai = atanf(fabs(vz / vx));
	float seita = atanf(fabs(vy / sqrtf(vx*vx+vz*vz)));
	*/
	float max_v;
	if(max_vx > max_vy)
		if(max_vx > max_vz)
			max_v = max_vx;
		else
			max_v = max_vz;
	else
		if(max_vy > max_vz)
			max_v = max_vy;
		else
			max_v = max_vz;
	if(max_v == 0)
		return;

	float dx = 0.5 * GRIDSIZE * vx/max_v;
	float dy = 0.5 * GRIDSIZE * vy/max_v;
	float dz = 0.5 * GRIDSIZE * vz/max_v;

	glBegin(GL_LINES);
	glColor3f(1.0f, 0.0f, 0.0f);
	glVertex3f(cx, cy, cz);
	glVertex3f(cx + dx, cy + dy, cz + dz);

	glColor3f(0.0f, 0.0f, 1.0f);
	glVertex3f(cx, cy, cz);
	glVertex3f(cx - dx, cy - dy, cz - dz);
	glEnd();
}

void FluidCube3D::draw_grid()
{
	glColor3f(0, 0, 0);
	glBegin(GL_LINES);
	for(int i = 0; i <= _N; i++)
		for(int j = 0; j <= _N; j++)
		{
			glVertex3f(i*GRIDSIZE, j*GRIDSIZE, 0);
			glVertex3f(i*GRIDSIZE, j*GRIDSIZE, LENGTH);
			glVertex3f(i*GRIDSIZE, 0, j*GRIDSIZE);
			glVertex3f(i*GRIDSIZE, LENGTH, j*GRIDSIZE);
			glVertex3f(0, i*GRIDSIZE, j*GRIDSIZE);
			glVertex3f(LENGTH, i*GRIDSIZE, j*GRIDSIZE);
		}
}

#endif