#include "stdafx.h"
#include "display.h"

int main(int argc, char* argv[])
{
	PRINT("Program Starting");
	glutInit(&argc, argv);
	//initialize the window 
	initialize();
	
	PRINT("Entering Main Loop");
	glutMainLoop(); //this starts the infinite loop
	PRINT("Exiting Program");

}
//1.  h2;  1/ dt
//2.  Store A