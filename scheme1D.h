/*
 * scheme1D.h
 *
 *  Created on: Mar 13, 2016
 *      Author: majun
 */

#ifndef SCHEME1D_H_
#define SCHEME1D_H_

#include "Matrix.h"
#include "Reconstruction.h"
#include <omp.h>
#include <cmath>
#include <iostream>

using namespace std;


class Grid1D
{
public:
	Grid1D(double IN_xInf,double IN_xSup,int IN_Ngrid);
	virtual ~Grid1D(void);
public:
	int Ngrid,Ncell;
	double xInf,xSup,dX,Lx;
	double* XG;
	double* XC;
};


class CSV1D
{
public:
	CSV1D(Grid1D& IN_GD);
	virtual ~CSV1D(void);

	void SetbyGrid(Grid1D& IN_GD);

	void Recon(string RecSch);
	void Flux_SplittingAdv(double Speed);
	void Diffusion(void);
	bool TimeAdvanceTEST(double timeRemain, double& timeInterval, double cfl,double Speed);
public:
	Grid1D* GD;

	double* Data;
	double* Data0;
	double* DtD;
	double* Fp;
	double* Fn;
	double* Flux;

	double* diff;

	int BoundaryCondition;
	double visc;
};



#endif /* SCHEME1D_H_ */
