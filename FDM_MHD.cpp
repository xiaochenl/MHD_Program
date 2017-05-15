//============================================================================
// Name        : FDM_MHD.cpp
// Author      : Trebu
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================


#include "Matrix.h"
#include <iostream>
#include "Task.h"
using namespace std;

int main() {
	cout << "!!!Hello World!!!" << endl; // prints !!!Hello World!!!

	//Test_Advection_WENO5();
	//TestAdvection2D();
	GLM_MHD2D_Run();
	return 0;
}
