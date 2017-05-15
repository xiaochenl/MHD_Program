/*
 * scheme1D.cpp
 *
 *  Created on: Mar 13, 2016
 *      Author: majun
 */

#include "scheme1D.h"

Grid1D::Grid1D(double IN_xInf,double IN_xSup,int IN_Ngrid)
{
	xInf = IN_xInf;
	xSup = IN_xSup;
	Ngrid = IN_Ngrid;

	Ncell = Ngrid-1;
	Lx = xSup-xInf;
	dX = Lx/Ncell;

	XG = new double [Ngrid];
	XC = new double [Ncell];

	for (int i=0; i<Ngrid; ++i) {
		XG[i] = xInf + i*dX;
	}
	for (int i=0; i<Ncell; ++i) {
		XC[i] = xInf + (i+0.5)*dX;
	}
}



Grid1D::~Grid1D(void)
{
	delete [] XG;
	delete [] XC;

	XG = NULL;
	XC = NULL;
}



CSV1D::CSV1D(Grid1D& IN_GD)
{
	GD = &IN_GD;

	Data = new double [IN_GD.Ngrid];
	Data0 = new double [IN_GD.Ngrid];
	DtD = new double [IN_GD.Ngrid];
	Fp = new double [IN_GD.Ngrid];
	Fn = new double [IN_GD.Ngrid];
	Flux = new double [IN_GD.Ngrid+1];
	diff = new double [IN_GD.Ngrid];

	BoundaryCondition = 0;
	visc = 0.0;
}


void CSV1D::SetbyGrid(Grid1D& IN_GD)
{
	GD = &IN_GD;

	Data = new double [IN_GD.Ngrid];
	Data0 = new double [IN_GD.Ngrid];
	DtD = new double [IN_GD.Ngrid];
	Fp = new double [IN_GD.Ngrid];
	Fn = new double [IN_GD.Ngrid];
	Flux = new double [IN_GD.Ngrid+1];
	diff = new double [IN_GD.Ngrid];

	BoundaryCondition = 0;
	visc = 0.0;
}



CSV1D::~CSV1D(void)
{
	delete [] Data;
	delete [] Data0;
	delete [] DtD;
	delete [] Fp;
	delete [] Fn;
	delete [] Flux;
	delete [] diff;

	Data = NULL;
	Data0 = NULL;
	DtD = NULL;
	Fp = NULL;
	Fn = NULL;
	Flux = NULL;
	diff = NULL;
}


void CSV1D::Recon(string RecSch)
{
	wenoKIND weno;
	if (RecSch=="WENO5p")	weno = WENO5;
	if (RecSch=="WENO5m")	weno = WENO5M;

#pragma omp parallel
	{
		TMP_WENO tmp;
#pragma omp for
		for (int i=3; i<(*GD).Ngrid-2; ++i) {
			Flux[i] = weno(Fp[i-3],Fp[i-2],Fp[i-1],Fp[i],Fp[i+1],tmp)
					+ weno(Fn[i+2],Fn[i+1],Fn[i],Fn[i-1],Fn[i-2],tmp);
		}
	}

	TMP_WENO tmp;
	if (0==BoundaryCondition) {			//	Periodical
		Flux[0] = weno(Fp[(*GD).Ngrid-4],Fp[(*GD).Ngrid-3],Fp[(*GD).Ngrid-2],Fp[0],Fp[1],tmp)
							+ weno(Fn[2],Fn[1],Fn[0],Fn[(*GD).Ngrid-2],Fn[(*GD).Ngrid-3],tmp);
		Flux[1] = weno(Fp[(*GD).Ngrid-3],Fp[(*GD).Ngrid-2],Fp[0],Fp[1],Fp[2],tmp)
							+ weno(Fn[2],Fn[1],Fn[0],Fn[(*GD).Ngrid-2],Fn[(*GD).Ngrid-3],tmp);
		Flux[2] = weno(Fp[(*GD).Ngrid-2],Fp[0],Fp[1],Fp[2],Fp[3],tmp)
							+ weno(Fn[4],Fn[3],Fn[2],Fn[1],Fn[0],tmp);
		Flux[(*GD).Ngrid-2] = weno(Fp[(*GD).Ngrid-5],Fp[(*GD).Ngrid-4],Fp[(*GD).Ngrid-3],Fp[(*GD).Ngrid-2],Fp[0],tmp)
							+ weno(Fn[1],Fn[0],Fn[(*GD).Ngrid-2],Fn[(*GD).Ngrid-3],Fn[(*GD).Ngrid-4],tmp);
		Flux[(*GD).Ngrid-1] = Flux[0];
		Flux[(*GD).Ngrid] = Flux[1];
	}
}


void CSV1D::Flux_SplittingAdv(double Speed)
{
	if (Speed>0.0) {
#pragma omp parallel for
		for (int _I = 0; _I < (*GD).Ngrid; ++_I) {
			Fp[_I] = Speed*Data[_I];
			Fn[_I] = 0.0;
		}
	} else {
#pragma omp parallel for
		for (int _I = 0; _I < (*GD).Ngrid; ++_I) {
			Fp[_I] = 0.0;
			Fn[_I] = Speed*Data[_I];
		}
	}
}



void CSV1D::Diffusion(void)
{
	int Ngrid = (*GD).Ngrid;
#pragma omp parallel for
	for (int _I = 3; _I < Ngrid-3; ++_I) {
		diff[_I] = -49.0/18.0*Data[_I] + 1.5*(Data[_I+1]+Data[_I-1]) - 0.15*(Data[_I+2]+Data[_I-2]) + 1.0/90.0*(Data[_I+3]+Data[_I-3]);
	}

	diff[2] = -49.0/18.0*Data[2] + 1.5*(Data[3]+Data[1]) - 0.15*(Data[4]+Data[0]) + 1.0/90.0*(Data[5]+Data[Ngrid-2]);
	diff[1] = -49.0/18.0*Data[1] + 1.5*(Data[2]+Data[0]) - 0.15*(Data[3]+Data[Ngrid-2]) + 1.0/90.0*(Data[4]+Data[Ngrid-3]);
	diff[0] = -49.0/18.0*Data[0] + 1.5*(Data[1]+Data[Ngrid-2]) - 0.15*(Data[2]+Data[Ngrid-3]) + 1.0/90.0*(Data[3]+Data[Ngrid-4]);
	diff[Ngrid-3] = -49.0/18.0*Data[Ngrid-3] + 1.5*(Data[Ngrid-2]+Data[Ngrid-4]) - 0.15*(Data[0]+Data[Ngrid-5]) + 1.0/90.0*(Data[1]+Data[Ngrid-6]);
	diff[Ngrid-2] = -49.0/18.0*Data[Ngrid-2] + 1.5*(Data[0]+Data[Ngrid-3]) - 0.15*(Data[1]+Data[Ngrid-4]) + 1.0/90.0*(Data[2]+Data[Ngrid-5]);
	diff[Ngrid-1] = diff[0];

#pragma omp parallel for
	for (int _I = 0; _I < Ngrid; ++_I) {
		diff[_I] /= (*GD).dX*(*GD).dX;
	}
}



bool CSV1D::TimeAdvanceTEST(double timeRemain, double& timeInterval, double cfl,double Speed)
{
	double a[3];
	a[0] = 0.0;
	a[1] = 0.75;
	a[2] = 1.0/3.0;
	bool reach = false;
	timeInterval = cfl*(*GD).dX/fabs(Speed);

	if (timeInterval - timeRemain >= 0.0) {
		timeInterval = timeRemain;
		reach = true;
	}

#pragma omp parallel for
	for (int _I = 0; _I < (*GD).Ngrid; ++_I) {
		Data0[_I] = Data[_I];
	}

	for (int step = 0; step < 3; ++step) {
		Diffusion();
		Flux_SplittingAdv(Speed);
		Recon("WENO5m");
#pragma omp parallel for
		for (int _I = 0; _I < (*GD).Ngrid; ++_I) {
			Data[_I] = a[step]*Data0[_I] + (1.0-a[step])*(Data[_I] + timeInterval*(-1.0/(*GD).dX*(Flux[_I+1]-Flux[_I]) + visc*diff[_I]));
		}
	}

	return reach;
}
