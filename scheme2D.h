/*
 * scheme2D.h
 *
 *  Created on: Mar 14, 2016
 *      Author: majun
 */
#pragma once
#ifndef SCHEME2D_H_
#define SCHEME2D_H_

#include "Matrix.h"
#include "Reconstruction.h"
#include <omp.h>
#include <cmath>
#include <iostream>

using namespace std;

class Grid2D
{
public:
	Grid2D(double IN_xInf,double IN_xSup,double IN_yInf,double IN_ySup,int IN_NGX,int IN_NGY);
	virtual ~Grid2D(void);
public:
	double xInf,xSup,yInf,ySup;
	int NGX,NGY,NCX,NCY;

	double dX,dY,Lx,Ly;

	double* XG;
	double* XC;
	double* YG;
	double* YC;
};

#endif /* SCHEME2D_H_ */


class CSV2D
{
public:
	CSV2D(void);
	CSV2D(Grid2D& IN_GD);
	virtual ~CSV2D(void) {}
	void SetByGrid(Grid2D& IN_GD);

	inline double* & operator [] (int row) const
	{
		return Data.m_RowStart[row];
	}

	void Recon(string RecSch);
	void FluxSplittingAdv(double Vx,double Vy);
	bool TimeAdvanceTEST(double timeRemain, double& timeInterval, double cfl,double Vx,double Vy);
public:
	Grid2D* GD;

	matrix<double> Data,Data0,DtD;
	matrix<double> FXp,FYp,FXn,FYn,FluxX,FluxY;

	int BoundaryCondition;
};



class GLM_MHD2D
{
public:
	GLM_MHD2D(Grid2D& IN_GD);
	virtual ~GLM_MHD2D(void);

	void Prim2CSV(void);
	void CSV2Prim(void);
	void CharaEst(void);
	void FluxSplitting(void);
	void GetDiffB(void);
	void SpacialDiscrt(string Recsch);

	bool TimeAdvance(double timeRemain, double& timeInterval, double cfl);
public:
	Grid2D* GD;
	CSV2D Mass,MomX,MomY,MomZ,Ener,Bx,By,Bz,Psi;
	CSV2D** CSV;
	matrix<double> Vx,Vy,Vz,Pres;
	matrix<double> Mass0,Vx0,Vy0,Vz0,Pres0,Bx0,By0,Bz0;
	matrix<double> LapBx,LapBy,LapBz;

	matrix<double> CharaX,CharaY;

	double gamma,gammah,eta,ch;
	int BoundaryCondition;

	int Nomp;
	double* maxVx;
	double* maxVy;
};

void GD2DParX(const matrix<double>& Sour,matrix<double>& Dest,double dX);
void GD2DParY(const matrix<double>& Sour,matrix<double>& Dest,double dY);
void GD2DParXX(const matrix<double>& Sour,matrix<double>& Dest,double dX);
void GD2DParYY(const matrix<double>& Sour,matrix<double>& Dest,double dY);
void GD2DLaplace(const matrix<double>& Sour,matrix<double>& Dest,double dX,double dY);
