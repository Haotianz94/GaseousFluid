#include "stdafx.h"
#include "display.h"

//#include <Eigen/Sparse>

int main(int argc, char* argv[])
{
	PRINT("Program Starting");

	glutInit(&argc, argv);
	//initialize the window 
	initialize();
	
	PRINT("Entering Main Loop");
	glutMainLoop(); //this starts the infinite loop
	PRINT("Exiting Program");
	return 0;

	/*Eigen::SparseMatrix<double, Eigen::RowMajor> mat(10, 10);         // default is column major
	mat.reserve(Eigen::VectorXi::Constant(10, 2));
	//mat.insert(i,j) = v_ij;                    // alternative: mat.coeffRef(i,j) += v_ij;
	//mat.makeCompressed();                        // optional
	std::cout<<mat<<std::endl;*/

	return 0;
}
//1.  h2;  1/ dt
//2.  Store A