#include "stdafx.h"
#include "fluidCube.h"
#include "display.h"
#include <math.h>
#include <ctime>
#include <cstdio>


FluidCube2D *cube;
int count = 0;
int wide, height;


void reshape(int _w, int _h){
	PRINT("start - reshape(" <<_w << "," << _h <<")");
	//glutPostRedisplay();
	wide = _w;
	height = _h;

	// Prevent a divide by zero, when window is too short
	// (you cant make a window of zero width).
	if (_h == 0)
		_h = 1;

	float ratio =  _w * 1.0 / _h;

	// Use the Projection Matrix
	glMatrixMode(GL_PROJECTION);

	// Reset Matrix
	glLoadIdentity();

	// Set the viewport to be the entire window
	glViewport(0, 0, _w, _h);

	// Set the correct perspective.
	gluPerspective(45,ratio,1,100);

	// Get Back to the Modelview
	glMatrixMode(GL_MODELVIEW);
	
	
	PRINT("done - reshape(" <<_w << "," << _h <<")");
}


int lastx, lasty;
clock_t start;

void mouseClick(int _button, int _state, int _x, int _y){
	PRINT("mouseClick(" <<_button << "," << _state << "," << _x << "," << _y <<")");
	if(_state)
	{
		PRINT("released");
		memset(cube->d0, 0, sizeof(float) * cube->size);
		memset(cube->Vx0, 0, sizeof(float) * cube->size);
		memset(cube->Vy0, 0, sizeof(float) * cube->size);

		//there may be some errors with the position
		int xx = _x / GRIDSIZE + 1;
		int yy = (height - _y) / GRIDSIZE ;
		clock_t time = clock() - start;
		float dxx = _x - lastx;
		float dyy = (height -_y) - lasty;
		float dxy = sqrtf(dxx*dxx + dyy*dyy);
		if(dxy == 0)
			return;

		cube->d0[IX(xx, yy)] = 10000;
		cube->Vx0[IX(xx, yy)] = dxx / dxy * time * DRAGSCALE;
		cube->Vy0[IX(xx, yy)] = dyy / dxy * time * DRAGSCALE;

		cube->simulate(false);
	}
	else
	{
		PRINT("clicked");
		lastx = _x;
		lasty = height -_y;
		start = clock();
	}
}


void keyEvent(unsigned char _key, int _x, int _y){
	
	PRINT("keyEvent(" << _key << "," << _x << "," << _y <<")");

	switch (_key)
	{
		case 27:  // ESC
			exit(0);
			break;
		case 'w':
			cube->setDisplayMode(VOTICITY);
			break;
		case 'l':
			cube->setDisplayMode(LIC);
			break;
		case 'n':
			cube->setDisplayMode(NONE);
			break;
		case 'd':
			cube->setDisplayMode(DENS);
			break;
		case 'v':
			cube->setDisplayMode(VELOCITY);
			break;
		default:
			break;
	}
}


void mouseDrag(int _x, int _y){
	PRINT("mouseDrag(" << _x << "," << _y <<"): displacement from click: (" << _x-lastx << "," << _y-lasty << ")");
}


void mouseMove(int _x, int _y){


}
void refresh(){

	PRINT("start - refresh()");
	cube->simulate(true);
	PRINT("done - refresh()");

}


void timer(int value) {

	if(count %FLOWTIME == 0)
	{
		memset(cube->d0, 0, sizeof(float) * cube->size);
		memset(cube->Vx0, 0, sizeof(float) * cube->size);
		memset(cube->Vy0, 0, sizeof(float) * cube->size);
		
		
		// for(int y = 1; y <= NUMGRIDH; y++)
		// {
		// 	cube->d0[IX(1, y)] = DENSITY;
		// 	cube->Vx0[IX(1, y)] = SPEED;  //10000~50000 for 2 vertexes
		// 	//cube->Vy0[IX(1, y)] = 0;
		// }
		
		float seita = (count % 30 + 75) / 180.0 * PI;
		int x = NUMGRIDW/2, y = 5;
		cube->d0[IX(x, y)] = DENSITY;
		cube->Vx0[IX(x, y)] = SPEED * cosf(seita);
		cube->Vy0[IX(x, y)] = SPEED * sinf(seita);
		
		cube->simulate(false);
	}
	else
	{
		glutPostRedisplay();
	}
	if(count < 111)//385 for vortex street 
		count ++;

	//glutPostRedisplay();
	glutTimerFunc(0, timer, 0); // next timer call milliseconds later
}


void initialize(){

#ifdef VISBLEW
	wide = VISBLEW;
#else
	wide = NUMGRIDW * GRIDSIZE + 20;
#endif

	height = NUMGRIDH * GRIDSIZE + 20;
	cube = new FluidCube2D(DIFFUSION, VISCOSITY, TIMESTEP);

	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE |GLUT_RGBA  | GLUT_STENCIL
                      | GLUT_ACCUM);
	glutInitWindowPosition(100, 100); //initial position of the window on the screen
	glutInitWindowSize(wide, height);
	glutCreateWindow("2D Fluid Simulation"); //create a window and set its focus

	//for the current window, set the callback functions:
	glutDisplayFunc(refresh); //infinite loop to draw on the window
	glutTimerFunc(0, timer, 0);
	//glutReshapeFunc(reshape); //called when the window is resized
	glutMouseFunc(mouseClick); //called when the mouse is clicked in the window
	glutKeyboardFunc(keyEvent); //called when a standard key is pressed
	//glutMotionFunc(mouseDrag); //called when the mouse is dragged after being clicked

	//glutIdleFunc(refresh);
	//glutPassiveMotionFunc(mouseMove); //called when the mouse moves (with/without click)
}