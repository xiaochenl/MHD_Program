/*
 * scheme2D.cpp
 *
 *  Created on: Mar 14, 2016
 *      Author: majun
 */


#include "scheme2D.h"


Grid2D::Grid2D(double IN_xInf,double IN_xSup,double IN_yInf,double IN_ySup,int IN_NGX,int IN_NGY)
{
	xInf = IN_xInf;
	xSup = IN_xSup;
	yInf = IN_yInf;
	ySup = IN_ySup;
	NGX = IN_NGX;
	NGY = IN_NGY;
	NCX = NGX-1;
	NCY = NGY-1;

	Lx = xSup-xInf;
	Ly = ySup-yInf;
	dX = Lx/NCX;
	dY = Ly/NCY;

	XG = new double [NGX];
	XC = new double [NCX];
	YG = new double [NGY];
	YC = new double [NCY];

	for (int i=0; i<NCX; i++) {
		XG[i] = xInf + i*dX;
		XC[i] = xInf + (i+0.5)*dX;
	}	XG[NCX] = xSup;
	for (int i=0; i<NCY; i++) {
		YG[i] = yInf + i*dY;
		YC[i] = yInf + (i+0.5)*dY;
	}	YG[NCY] = ySup;
}


Grid2D::~Grid2D(void)
{
	delete [] XG;
	delete [] XC;
	delete [] YG;
	delete [] YC;

	XG = NULL;
	XC = NULL;
	YG = NULL;
	YC = NULL;
}



CSV2D::CSV2D(void)
{
	GD = NULL;
	BoundaryCondition = 0;
}



CSV2D::CSV2D(Grid2D& IN_GD)
{
	GD = &IN_GD;

	Data.Resize(IN_GD.NGY,IN_GD.NGX);
	Data0.Resize(IN_GD.NGY,IN_GD.NGX);
	DtD.Resize(IN_GD.NGY,IN_GD.NGX);
	FXp.Resize(IN_GD.NGY,IN_GD.NGX);
	FXn.Resize(IN_GD.NGY,IN_GD.NGX);
	FYp.Resize(IN_GD.NGY,IN_GD.NGX);
	FYn.Resize(IN_GD.NGY,IN_GD.NGX);
	FluxX.Resize(IN_GD.NGY,IN_GD.NGX+1);
	FluxY.Resize(IN_GD.NGY+1,IN_GD.NGX);

	BoundaryCondition = 0;
}


void CSV2D::SetByGrid(Grid2D& IN_GD)
{
	GD = &IN_GD;

	Data.Resize(IN_GD.NGY,IN_GD.NGX);
	Data0.Resize(IN_GD.NGY,IN_GD.NGX);
	DtD.Resize(IN_GD.NGY,IN_GD.NGX);
	FXp.Resize(IN_GD.NGY,IN_GD.NGX);
	FXn.Resize(IN_GD.NGY,IN_GD.NGX);
	FYp.Resize(IN_GD.NGY,IN_GD.NGX);
	FYn.Resize(IN_GD.NGY,IN_GD.NGX);
	FluxX.Resize(IN_GD.NGY,IN_GD.NGX+1);
	FluxY.Resize(IN_GD.NGY+1,IN_GD.NGX);

	BoundaryCondition = 0;
}


void CSV2D::Recon(string RecSch)
{
	int Nx = (*GD).NGX;
	int Ny = (*GD).NGY;

	wenoKIND weno;
	if (RecSch=="WENO5p")	weno = WENO5;
	if (RecSch=="WENO5m")	weno = WENO5M;

#pragma omp parallel
	{
		TMP_WENO tmp;
#pragma omp for
		for (int j=0; j<Ny; ++j) {
			for (int i=3; i<Nx-2; ++i) {
				FluxX[j][i] = weno(FXp[j][i-3],FXp[j][i-2],FXp[j][i-1],FXp[j][i],FXp[j][i+1],tmp)
						+ weno(FXn[j][i+2],FXn[j][i+1],FXn[j][i],FXn[j][i-1],FXn[j][i-2],tmp);
			}
			FluxX[j][0] = weno(FXp[j][Nx-4],FXp[j][Nx-3],FXp[j][Nx-2],FXp[j][0],FXp[j][1],tmp)
								+ weno(FXn[j][2],FXn[j][1],FXn[j][0],FXn[j][Nx-2],FXn[j][Nx-3],tmp);
			FluxX[j][1] = weno(FXp[j][Nx-3],FXp[j][Nx-2],FXp[j][0],FXp[j][1],FXp[j][2],tmp)
								+ weno(FXn[j][3],FXn[j][2],FXn[j][1],FXn[j][0],FXn[j][Nx-2],tmp);
			FluxX[j][2] = weno(FXp[j][Nx-2],FXp[j][0],FXp[j][1],FXp[j][2],FXp[j][3],tmp)
								+ weno(FXn[j][4],FXn[j][3],FXn[j][2],FXn[j][1],FXn[j][0],tmp);
			FluxX[j][Nx-2] = weno(FXp[j][Nx-5],FXp[j][Nx-4],FXp[j][Nx-3],FXp[j][Nx-2],FXp[j][0],tmp)
								+ weno(FXn[j][1],FXn[j][0],FXn[j][Nx-2],FXn[j][Nx-3],FXn[j][Nx-4],tmp);
			FluxX[j][Nx-1] = FluxX[j][0];
			FluxX[j][Nx] = FluxX[j][1];
		}
		//
#pragma omp for
		for (int i=0; i<Nx; ++i) {
			for (int j=3; j<Ny-2; ++j) {
				FluxY[j][i] = weno(FYp[j-3][i],FYp[j-2][i],FYp[j-1][i],FYp[j][i],FYp[j+1][i],tmp)
						+ weno(FYn[j+2][i],FYn[j+1][i],FYn[j][i],FYn[j-1][i],FYn[j-2][i],tmp);
			}
			FluxY[0][i] = weno(FYp[Ny-4][i],FYp[Ny-3][i],FYp[Ny-2][i],FYp[0][i],FYp[1][i],tmp)
								+ weno(FYn[2][i],FYn[1][i],FYn[0][i],FYn[Ny-2][i],FYn[Ny-3][i],tmp);
			FluxY[1][i] = weno(FYp[Ny-3][i],FYp[Ny-2][i],FYp[0][i],FYp[1][i],FYp[2][i],tmp)
								+ weno(FYn[3][i],FYn[2][i],FYn[1][i],FYn[0][i],FYn[Ny-2][i],tmp);
			FluxY[2][i] = weno(FYp[Ny-2][i],FYp[0][i],FYp[1][i],FYp[2][i],FYp[3][i],tmp)
								+ weno(FYn[4][i],FYn[3][i],FYn[2][i],FYn[1][i],FYn[0][i],tmp);
			FluxY[Ny-2][i] = weno(FYp[Ny-5][i],FYp[Ny-4][i],FYp[Ny-3][i],FYp[Ny-2][i],FYp[0][i],tmp)
								+ weno(FYn[1][i],FYn[0][i],FYn[Ny-2][i],FYn[Ny-3][i],FYn[Ny-4][i],tmp);
			FluxY[Ny-1][i] = FluxY[0][i];
			FluxY[Ny][i] = FluxY[1][i];
		}
	}
}


void CSV2D::FluxSplittingAdv(double Vx,double Vy)
{
	int Nx = (*GD).NGX;
	int Ny = (*GD).NGY;
	if (Vx>0.0) {
#pragma omp parallel for
		for (int j=0; j<Ny; ++j) {
			for (int i=0; i<Nx; ++i) {
				FXp[j][i] = Vx*Data[j][i];
				FXn[j][i] = 0.0;
			}
		}
	} else {
#pragma omp parallel for
		for (int j=0; j<Ny; ++j) {
			for (int i=0; i<Nx; ++i) {
				FXp[j][i] = 0.0;
				FXn[j][i] = Vx*Data[j][i];
			}
		}
	}

	if (Vy>0.0) {
#pragma omp parallel for
		for (int j=0; j<Ny; ++j) {
			for (int i=0; i<Nx; ++i) {
				FYp[j][i] = Vy*Data[j][i];
				FYn[j][i] = 0.0;
			}
		}
	} else {
#pragma omp parallel for
		for (int j=0; j<Ny; ++j) {
			for (int i=0; i<Nx; ++i) {
				FYp[j][i] = 0.0;
				FYn[j][i] = Vy*Data[j][i];
			}
		}
	}
}


bool CSV2D::TimeAdvanceTEST(double timeRemain, double& timeInterval, double cfl,double Vx,double Vy)
{
	bool reach = false;
	double a[3];
	a[0] = 0.0;
	a[1] = 0.75;
	a[2] = 1.0/3.0;

	timeInterval = 0.5*cfl*min((*GD).dY/fabs(Vy),(*GD).dX/fabs(Vx));
	timeInterval = min(timeRemain,timeInterval);
	if (timeInterval - timeRemain >= 0.0) {
		timeInterval = timeRemain;
		reach = true;
	}

#pragma omp parallel for
	for (int i=0; i<(*GD).NGX; ++i) {
		for (int j=0; j<(*GD).NGY; ++j) {
			Data0[j][i] = Data[j][i];
		}
	}

	for (int step = 0; step < 3; ++step) {
		//Diffusion();
		FluxSplittingAdv(Vx,Vy);
		Recon("WENO5m");
#pragma omp parallel for
		for (int i=0; i<(*GD).NGX; ++i) {
			for (int j=0; j<(*GD).NGY; ++j) {
				Data[j][i] = a[step]*Data0[j][i] + (1.0-a[step])*(Data[j][i] + timeInterval*(-1.0/(*GD).dX*(FluxX[j][i+1]-FluxX[j][i])-1.0/(*GD).dY*(FluxY[j+1][i]-FluxY[j][i]) ) );
			}
		}
	}

	return reach;
}




GLM_MHD2D::GLM_MHD2D(Grid2D& IN_GD)
{
	GD = &IN_GD;
	Mass.SetByGrid(IN_GD);
	MomX.SetByGrid(IN_GD);
	MomY.SetByGrid(IN_GD);
	MomZ.SetByGrid(IN_GD);
	Ener.SetByGrid(IN_GD);
	Bx.SetByGrid(IN_GD);
	By.SetByGrid(IN_GD);
	Bz.SetByGrid(IN_GD);
	Psi.SetByGrid(IN_GD);

	CSV = new CSV2D* [9];
	CSV[0] = &Mass;
	CSV[1] = &MomX;
	CSV[2] = &MomY;
	CSV[3] = &MomZ;
	CSV[4] = &Ener;
	CSV[5] = &Bx;
	CSV[6] = &By;
	CSV[7] = &Bz;
	CSV[8] = &Psi;

	Vx.Resize(IN_GD.NGY,IN_GD.NGX);
	Vy.Resize(IN_GD.NGY,IN_GD.NGX);
	Vz.Resize(IN_GD.NGY,IN_GD.NGX);
	Pres.Resize(IN_GD.NGY,IN_GD.NGX);

	Mass0.Resize(IN_GD.NGY,IN_GD.NGX);
	Vx0.Resize(IN_GD.NGY,IN_GD.NGX);
	Vy0.Resize(IN_GD.NGY,IN_GD.NGX);
	Vz0.Resize(IN_GD.NGY,IN_GD.NGX);
	Pres0.Resize(IN_GD.NGY,IN_GD.NGX);
	Bx0.Resize(IN_GD.NGY,IN_GD.NGX);
	By0.Resize(IN_GD.NGY,IN_GD.NGX);
	Bz0.Resize(IN_GD.NGY,IN_GD.NGX);

	LapBx.Resize(IN_GD.NGY,IN_GD.NGX);
	LapBy.Resize(IN_GD.NGY,IN_GD.NGX);
	LapBz.Resize(IN_GD.NGY,IN_GD.NGX);

	CharaX.Resize(IN_GD.NGY,IN_GD.NGX);
	CharaY.Resize(IN_GD.NGY,IN_GD.NGX);

	gamma = 5.0/3.0;
	gammah = gamma-1.0;
	eta = 0.0;
	ch = 0.0;
	BoundaryCondition = 0;

#pragma omp parallel
	{
		if (omp_get_thread_num()==0)	Nomp = omp_get_num_threads();
	}

	maxVx = new double [Nomp];
	maxVy = new double [Nomp];
}


GLM_MHD2D::~GLM_MHD2D(void)
{
	delete [] CSV;
	delete [] maxVx;
	delete [] maxVy;
	CSV = NULL;
	maxVx = NULL;
	maxVy = NULL;
}


void GLM_MHD2D::Prim2CSV(void)
{
	int Nx = (*GD).NGX;
	int Ny = (*GD).NGY;
#pragma omp parallel for
	for (int i=0; i<Ny; ++i) {
		for (int j=0; j<Nx; ++j) {
			MomX.Data[i][j] = Mass.Data[i][j]*Vx0[i][j] + Mass0[i][j]*Vx[i][j] + Mass.Data[i][j]*Vx[i][j];
			MomY.Data[i][j] = Mass.Data[i][j]*Vy0[i][j] + Mass0[i][j]*Vy[i][j] + Mass.Data[i][j]*Vy[i][j];
			MomZ.Data[i][j] = Mass.Data[i][j]*Vz0[i][j] + Mass0[i][j]*Vz[i][j] + Mass.Data[i][j]*Vz[i][j];
//			Ener.Data[i][j] = Pres[i][j]/gammah + 0.5*Mass.Data[i][j]*(Vx0[i][j]*Vx0[i][j]+Vy0[i][j]*Vy0[i][j]+Vz0[i][j]*Vz0[i][j])
//					        + (Mass0[i][j]+Mass.Data[i][j])*( (Vx0[i][j]*Vx[i][j]+Vy0[i][j]*Vy[i][j]+Vz0[i][j]*Vz[i][j])
//					        + 0.5*(Vx[i][j]*Vx[i][j]+Vy[i][j]*Vy[i][j]+Vz[i][j]*Vz[i][j]) )
//							+ Bx0[i][j]*Bx.Data[i][j]+By0[i][j]*By.Data[i][j]+Bz0[i][j]*Bz.Data[i][j]
//							+ 0.5*(Bx.Data[i][j]*Bx.Data[i][j]+By.Data[i][j]*By.Data[i][j]+Bz.Data[i][j]*Bz.Data[i][j]);
			Ener.Data[i][j] = Pres[i][j] + gammah*(0.5*Mass.Data[i][j]*(Vx0[i][j]*Vx0[i][j]+Vy0[i][j]*Vy0[i][j]+Vz0[i][j]*Vz0[i][j])
					        + (Mass0[i][j]+Mass.Data[i][j])*( (Vx0[i][j]*Vx[i][j]+Vy0[i][j]*Vy[i][j]+Vz0[i][j]*Vz[i][j])
					        + 0.5*(Vx[i][j]*Vx[i][j]+Vy[i][j]*Vy[i][j]+Vz[i][j]*Vz[i][j]) )
							+ Bx0[i][j]*Bx.Data[i][j]+By0[i][j]*By.Data[i][j]+Bz0[i][j]*Bz.Data[i][j]
							+ 0.5*(Bx.Data[i][j]*Bx.Data[i][j]+By.Data[i][j]*By.Data[i][j]+Bz.Data[i][j]*Bz.Data[i][j]));

		}
	}
}


void GLM_MHD2D::CSV2Prim(void)
{
	int Nx = (*GD).NGX;
	int Ny = (*GD).NGY;
#pragma omp parallel for
	for (int i=0; i<Ny; ++i) {
		for (int j=0; j<Nx; ++j) {
			Vx[i][j] = (MomX.Data[i][j]-Mass.Data[i][j]*Vx0[i][j])/(Mass.Data[i][j]+Mass0[i][j]);
			Vy[i][j] = (MomY.Data[i][j]-Mass.Data[i][j]*Vy0[i][j])/(Mass.Data[i][j]+Mass0[i][j]);
			Vz[i][j] = (MomZ.Data[i][j]-Mass.Data[i][j]*Vz0[i][j])/(Mass.Data[i][j]+Mass0[i][j]);
			Pres[i][j] = Ener.Data[i][j]-gammah*(0.5*Mass.Data[i][j]*(Vx0[i][j]*Vx0[i][j]+Vy0[i][j]*Vy0[i][j]+Vz0[i][j]*Vz0[i][j])
			        		+ (Mass0[i][j]+Mass.Data[i][j])*( (Vx0[i][j]*Vx[i][j]+Vy0[i][j]*Vy[i][j]+Vz0[i][j]*Vz[i][j])
			        		+ 0.5*(Vx[i][j]*Vx[i][j]+Vy[i][j]*Vy[i][j]+Vz[i][j]*Vz[i][j]) )
							+ Bx0[i][j]*Bx.Data[i][j]+By0[i][j]*By.Data[i][j]+Bz0[i][j]*Bz.Data[i][j]
							+ 0.5*(Bx.Data[i][j]*Bx.Data[i][j]+By.Data[i][j]*By.Data[i][j]+Bz.Data[i][j]*Bz.Data[i][j]));
		}
	}
}


void GLM_MHD2D::CharaEst(void)
{
	int Nx = (*GD).NGX;
	int Ny = (*GD).NGY;
	double pp,pb,px,py,pz;
#pragma omp parallel for private(pp,pb,px,py,pz)
	for (int i=0; i<Ny; ++i) {
		for (int j=0; j<Nx; ++j) {
			pp = gamma*(Pres0[i][j]+Pres[i][j]);
			px = (Bx.Data[i][j]+Bx0[i][j])*(Bx.Data[i][j]+Bx0[i][j]);
			py = (By.Data[i][j]+By0[i][j])*(By.Data[i][j]+By0[i][j]);
			pz = (Bz.Data[i][j]+Bz0[i][j])*(Bz.Data[i][j]+Bz0[i][j]);
			pb = px+py+pz;

			CharaX[i][j] = sqrt(0.5*(pp+pb+sqrt((pp+pb)*(pp+pb)-4*pp*px))/(Mass.Data[i][j]+Mass0[i][j]))+fabs(Vx[i][j]+Vx0[i][j]);
			CharaY[i][j] = sqrt(0.5*(pp+pb+sqrt((pp+pb)*(pp+pb)-4*pp*py))/(Mass.Data[i][j]+Mass0[i][j]))+fabs(Vy[i][j]+Vy0[i][j]);
		}
	}

#pragma omp parallel
	{
		int id = omp_get_thread_num();
		maxVx[id] = 0.0;
		maxVy[id] = 0.0;
		for (int i=id; i<Ny; i+=Nomp) {
			for (int j=0; j<Nx; ++j) {
				if (maxVx[id]<CharaX[i][j])	maxVx[id] = CharaX[i][j];
				if (maxVy[id]<CharaY[i][j])	maxVy[id] = CharaY[i][j];
			}
		}
	}

	ch = 0.0;
	for (int i=0; i<Nomp; ++i) {
		if (ch<maxVx[i])	ch = maxVx[i];
		if (ch<maxVy[i])	ch = maxVy[i];
	}

}


void GLM_MHD2D::FluxSplitting(void)
{
	int Nx = (*GD).NGX;
	int Ny = (*GD).NGY;
	double H0,H1,mu0,mu1;
#pragma omp parallel for private(H0,H1,mu0,mu1)
	for (int i=0; i<Ny; ++i) {
		for (int j=0; j<Nx; ++j) {
			H0 = gamma*Pres0[i][j]+gammah*(0.5*Mass0[i][j]*(Vx0[i][j]*Vx0[i][j]+Vy0[i][j]*Vy0[i][j]+Vz0[i][j]*Vz0[i][j])
					         + Bx0[i][j]*Bx0[i][j]+By0[i][j]*By0[i][j]+Bz0[i][j]*Bz0[i][j]);
			H1 = Ener.Data[i][j] + gammah*(Pres[i][j] + Bx0[i][j]*Bx.Data[i][j] + By0[i][j]*By.Data[i][j] + Bz0[i][j]*Bz.Data[i][j]
							 + 0.5*(Bx.Data[i][j]*Bx.Data[i][j]+By.Data[i][j]*By.Data[i][j]+Bz.Data[i][j]*Bz.Data[i][j]) );
			mu0 = gammah*(Bx0[i][j]*Vx0[i][j]+By0[i][j]*Vy0[i][j]+Bz0[i][j]*Vz0[i][j]);
			mu1 = gammah*(Bx0[i][j]*Vx[i][j]+By0[i][j]*Vy[i][j]+Bz0[i][j]*Vz[i][j]  + Bx.Data[i][j]*Vx[i][j]+By.Data[i][j]*Vy[i][j]+Bz.Data[i][j]*Vz[i][j]
				+ Bx.Data[i][j]*Vx0[i][j]+By.Data[i][j]*Vy0[i][j]+Bz.Data[i][j]*Vz0[i][j]);
			///////////////////////////////////////////////////////////////////////////////
			Mass.FluxX[i][j] = MomX.Data[i][j];
			MomX.FluxX[i][j] = Pres[i][j] + Mass.Data[i][j]*Vx0[i][j]*Vx0[i][j]
					         + (Mass0[i][j]+Mass.Data[i][j])*Vx[i][j]*(Vx[i][j]+2.0*Vx0[i][j])
							 + By.Data[i][j]*By0[i][j] + Bz.Data[i][j]*Bz0[i][j] - Bx.Data[i][j]*Bx0[i][j]
							 + 0.5*(-Bx.Data[i][j]*Bx.Data[i][j]+By.Data[i][j]*By.Data[i][j]+Bz.Data[i][j]*Bz.Data[i][j]);
			MomY.FluxX[i][j] = Mass0[i][j]*Vx0[i][j]*Vy[i][j]+Mass0[i][j]*Vx[i][j]*Vy0[i][j]+Mass.Data[i][j]*Vx0[i][j]*Vy0[i][j]
							 + Mass0[i][j]*Vx[i][j]*Vy[i][j]+Mass.Data[i][j]*Vx0[i][j]*Vy[i][j]+Mass.Data[i][j]*Vx[i][j]*Vy0[i][j]
							 + Mass.Data[i][j]*Vx[i][j]*Vy[i][j]
							 - Bx0[i][j]*By.Data[i][j]-Bx.Data[i][j]*By0[i][j]-Bx.Data[i][j]*By.Data[i][j];
			MomZ.FluxX[i][j] = Mass0[i][j]*Vx0[i][j]*Vz[i][j]+Mass0[i][j]*Vx[i][j]*Vz0[i][j]+Mass.Data[i][j]*Vx0[i][j]*Vz0[i][j]
							 + Mass0[i][j]*Vx[i][j]*Vz[i][j]+Mass.Data[i][j]*Vx0[i][j]*Vz[i][j]+Mass.Data[i][j]*Vx[i][j]*Vz0[i][j]
							 + Mass.Data[i][j]*Vx[i][j]*Vz[i][j]
							 - Bx0[i][j]*Bz.Data[i][j]-Bx.Data[i][j]*Bz0[i][j]-Bx.Data[i][j]*Bz.Data[i][j];
			Ener.FluxX[i][j] = Vx[i][j]*H0+Vx0[i][j]*H1+Vx[i][j]*H1-Bx.Data[i][j]*mu0-Bx0[i][j]*mu1-Bx.Data[i][j]*mu1;
			Bx.FluxX[i][j] = Psi.Data[i][j];
			By.FluxX[i][j] = By0[i][j]*Vx[i][j]+By.Data[i][j]*Vx0[i][j]-Bx0[i][j]*Vy[i][j]
						   - Bx.Data[i][j]*Vy0[i][j]+By.Data[i][j]*Vx[i][j]-Bx.Data[i][j]*Vy[i][j];
			Bz.FluxX[i][j] = Bz0[i][j]*Vx[i][j]+Bz.Data[i][j]*Vx0[i][j]-Bx0[i][j]*Vz[i][j]
						   - Bx.Data[i][j]*Vz0[i][j]+Bz.Data[i][j]*Vx[i][j]-Bx.Data[i][j]*Vz[i][j];
			Psi.FluxX[i][j] = ch*ch*Bx.Data[i][j];
			///////////////////////////////////////////////////////////////////////////////
			Mass.FluxY[i][j] = MomY.Data[i][j];
			MomX.FluxY[i][j] = MomY.FluxX[i][j];
			MomY.FluxY[i][j] = Pres[i][j] + Mass.Data[i][j]*Vy0[i][j]*Vy0[i][j]
							 + (Mass0[i][j]+Mass.Data[i][j])*Vy[i][j]*(Vy[i][j]+2.0*Vy0[i][j])
							 + Bx.Data[i][j]*Bx0[i][j] + Bz.Data[i][j]*Bz0[i][j] - By.Data[i][j]*By0[i][j]
							 + 0.5*(-By.Data[i][j]*By.Data[i][j]+Bx.Data[i][j]*Bx.Data[i][j]+Bz.Data[i][j]*Bz.Data[i][j]);
			MomZ.FluxY[i][j] = Mass0[i][j]*Vy0[i][j]*Vz[i][j]+Mass0[i][j]*Vy[i][j]*Vz0[i][j]+Mass.Data[i][j]*Vy0[i][j]*Vz0[i][j]
							 + Mass0[i][j]*Vy[i][j]*Vz[i][j]+Mass.Data[i][j]*Vy0[i][j]*Vz[i][j]+Mass.Data[i][j]*Vy[i][j]*Vz0[i][j]
							 + Mass.Data[i][j]*Vz[i][j]*Vy[i][j]
							 - By0[i][j]*Bz.Data[i][j]-By.Data[i][j]*Bz0[i][j]-By.Data[i][j]*Bz.Data[i][j];
			Ener.FluxY[i][j] = Vy[i][j]*H0+Vy0[i][j]*H1+Vy[i][j]*H1-By.Data[i][j]*mu0-By0[i][j]*mu1-By.Data[i][j]*mu1;
			Bx.FluxY[i][j] = -By.FluxX[i][j];
			By.FluxY[i][j] = Psi.Data[i][j];
			Bz.FluxY[i][j] = Bz0[i][j]*Vy[i][j]+Bz.Data[i][j]*Vy0[i][j]-By0[i][j]*Vz[i][j]
							 - By.Data[i][j]*Vz0[i][j]+Bz.Data[i][j]*Vy[i][j]-By.Data[i][j]*Vz[i][j];
			Psi.FluxY[i][j] = ch*ch*By.Data[i][j];
			///////////////////////////////////////////////////////////////////////////////
			for (int k=0; k<9; ++k) {
				(*CSV[k]).FXp[i][j] = 0.5*((*CSV[k]).FluxX[i][j]+CharaX[i][j]*(*CSV[k]).Data[i][j]);
				(*CSV[k]).FXn[i][j] = 0.5*((*CSV[k]).FluxX[i][j]-CharaX[i][j]*(*CSV[k]).Data[i][j]);
				(*CSV[k]).FYp[i][j] = 0.5*((*CSV[k]).FluxY[i][j]+CharaY[i][j]*(*CSV[k]).Data[i][j]);
				(*CSV[k]).FYn[i][j] = 0.5*((*CSV[k]).FluxY[i][j]-CharaY[i][j]*(*CSV[k]).Data[i][j]);
//				(*CSV[k]).FXp[i][j] = 0.5*((*CSV[k]).FluxX[i][j]+ch*(*CSV[k]).Data[i][j]);
//				(*CSV[k]).FXn[i][j] = 0.5*((*CSV[k]).FluxX[i][j]-ch*(*CSV[k]).Data[i][j]);
//				(*CSV[k]).FYp[i][j] = 0.5*((*CSV[k]).FluxY[i][j]+ch*(*CSV[k]).Data[i][j]);
//				(*CSV[k]).FYn[i][j] = 0.5*((*CSV[k]).FluxY[i][j]-ch*(*CSV[k]).Data[i][j]);
				//(*CSV[k]).FluxX[i][j] = 0.0;
				//(*CSV[k]).FluxY[i][j] = 0.0;
			}
		}
	}
}


void GLM_MHD2D::SpacialDiscrt(string Recsch)
{
	CharaEst();
	FluxSplitting();
	for (int k=0; k<9; ++k) {
		(*CSV[k]).Recon(Recsch);
	}
}


bool GLM_MHD2D::TimeAdvance(double timeRemain, double& timeInterval, double cfl)
{
	double etaThresh = 1.0e-16;
	bool reach = false;
	double a[3];
	a[0] = 0.0;
	a[1] = 0.75;
	a[2] = 1.0/3.0;
	string Recsch = "WENO5m";


	SpacialDiscrt(Recsch);
	timeInterval = 0.5*cfl/ch*min((*GD).dX,(*GD).dY);

	if (timeInterval-timeRemain>=0.0) {
		timeInterval = timeRemain;
		reach = true;
	}

	double decayRate = exp(-timeInterval*ch/0.18);

#pragma omp parallel for
	for (int i=0; i<(*GD).NGY; ++i) {
		for (int j=0; j<(*GD).NGX; ++j)	{
			for (int k=0; k<9; ++k) {
				(*CSV[k]).Data0[i][j] = (*CSV[k]).Data[i][j];
			}
		}
	}

#pragma omp parallel for
	for (int i=0; i<(*GD).NGY; ++i) {
		for (int j=0; j<(*GD).NGX; ++j)	{
			for (int k=0; k<9; ++k) {
				(*CSV[k]).Data[i][j] += timeInterval*(-1.0/(*GD).dX*((*CSV[k]).FluxX[i][j+1]
									   -(*CSV[k]).FluxX[i][j])-1.0/(*GD).dY*((*CSV[k]).FluxY[i+1][j]-(*CSV[k]).FluxY[i][j]) );
			}
		}
	}

	//	Currently valid for gamma = 1.0 and not too small eta
	if (eta > etaThresh) {
		GetDiffB();
#pragma omp parallel for
		for (int i=0; i<(*GD).NGY; ++i) {
			for (int j=0; j<(*GD).NGX; ++j)	{
				Bx.Data[i][j] += timeInterval*eta*LapBx[i][j];
				By.Data[i][j] += timeInterval*eta*LapBy[i][j];
				Bz.Data[i][j] += timeInterval*eta*LapBz[i][j];
			}
		}
	}
	CSV2Prim();
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	for (int step = 1; step < 3; ++step) {
		SpacialDiscrt(Recsch);
#pragma omp parallel for
		for (int i=0; i<(*GD).NGY; ++i) {
			for (int j=0; j<(*GD).NGX; ++j) {
				for (int k=0; k<9; ++k) {
					(*CSV[k]).Data[i][j] = a[step]*(*CSV[k]).Data0[i][j] + (1.0-a[step])*((*CSV[k]).Data[i][j]
										 + timeInterval*(-1.0/(*GD).dX*((*CSV[k]).FluxX[i][j+1]-(*CSV[k]).FluxX[i][j])
										 - 1.0/(*GD).dY*((*CSV[k]).FluxY[i+1][j]-(*CSV[k]).FluxY[i][j]) ) );
				}
			}
		}
		//
		//	Currently valid for gamma = 1.0 and not too small eta
		if (eta > etaThresh) {
			GetDiffB();
#pragma omp parallel for
			for (int i=0; i<(*GD).NGY; ++i) {
				for (int j=0; j<(*GD).NGX; ++j)	{
					Bx.Data[i][j] += (1.0-a[step])*timeInterval*eta*LapBx[i][j];
					By.Data[i][j] += (1.0-a[step])*timeInterval*eta*LapBy[i][j];
					Bz.Data[i][j] += (1.0-a[step])*timeInterval*eta*LapBz[i][j];
				}
			}
		}

		CSV2Prim();
	}

#pragma omp parallel for
	for (int i=0; i<(*GD).NGY; ++i) {
		for (int j=0; j<(*GD).NGX; ++j) {
			Psi.Data[i][j] *= decayRate;
		}
	}

	return reach;
}


void GLM_MHD2D::GetDiffB(void)
{
	GD2DLaplace(Bx.Data,LapBx,(*GD).dX,(*GD).dY);
	GD2DLaplace(By.Data,LapBy,(*GD).dX,(*GD).dY);
	GD2DLaplace(Bz.Data,LapBz,(*GD).dX,(*GD).dY);
}



void GD2DParX(const matrix<double>& Sour,matrix<double>& Dest,double dX)
{
	int Nx = Sour.m_Column;
	int Ny = Sour.m_Row;
#pragma omp parallel for
	for (int i=0; i<Ny; ++i) {
		for (int j=3; j<Nx-3; ++j) {
			Dest[i][j] = 0.75*(Sour[i][j+1]-Sour[i][j-1]) - 0.15*(Sour[i][j+2]-Sour[i][j-2]) + 1.0/60.0*(Sour[i][j+3]-Sour[i][j-3]);
		}
		Dest[i][2] = 0.75*(Sour[i][3]-Sour[i][1]) - 0.15*(Sour[i][4]-Sour[i][0]) + 1.0/60.0*(Sour[i][5]-Sour[i][Nx-2]);
		Dest[i][1] = 0.75*(Sour[i][2]-Sour[i][0]) - 0.15*(Sour[i][3]-Sour[i][Nx-2]) + 1.0/60.0*(Sour[i][4]-Sour[i][Nx-3]);
		Dest[i][0] = 0.75*(Sour[i][1]-Sour[i][Nx-2]) - 0.15*(Sour[i][2]-Sour[i][Nx-3]) + 1.0/60.0*(Sour[i][3]-Sour[i][Nx-4]);
		Dest[i][Nx-3] = 0.75*(Sour[i][Nx-2]-Sour[i][Nx-4]) - 0.15*(Sour[i][0]-Sour[i][Nx-5]) + 1.0/60.0*(Sour[i][1]-Sour[i][Nx-6]);
		Dest[i][Nx-2] = 0.75*(Sour[i][0]-Sour[i][Nx-3]) - 0.15*(Sour[i][1]-Sour[i][Nx-4]) + 1.0/60.0*(Sour[i][2]-Sour[i][Nx-5]);
		Dest[i][Nx-1] = Dest[i][0];
		//
		for (int j=0; j<Nx; ++j)		Dest[i][j] *= 1.0/dX;
	}
}



void GD2DParY(const matrix<double>& Sour,matrix<double>& Dest,double dY)
{
	int Nx = Sour.m_Column;
	int Ny = Sour.m_Row;
#pragma omp parallel for
	for (int j=0; j<Nx; ++j) {
		for (int i=3; i<Ny-3; ++i) {
			Dest[i][j] = 0.75*(Sour[i+1][j]-Sour[i-1][j]) - 0.15*(Sour[i+2][j]-Sour[i-2][j]) + 1.0/60.0*(Sour[i+3][j]-Sour[i-3][j]);
		}
		Dest[2][j] = 0.75*(Sour[3][j]-Sour[1][j]) - 0.15*(Sour[4][j]-Sour[0][j]) + 1.0/60.0*(Sour[5][j]-Sour[Ny-2][j]);
		Dest[1][j] = 0.75*(Sour[2][j]-Sour[0][j]) - 0.15*(Sour[3][j]-Sour[Ny-2][j]) + 1.0/60.0*(Sour[4][j]-Sour[Ny-3][j]);
		Dest[0][j] = 0.75*(Sour[1][j]-Sour[Ny-2][j]) - 0.15*(Sour[2][j]-Sour[Ny-3][j]) + 1.0/60.0*(Sour[3][j]-Sour[Ny-4][j]);
		Dest[Ny-3][j] = 0.75*(Sour[Ny-2][j]-Sour[Ny-4][j]) - 0.15*(Sour[0][j]-Sour[Ny-5][j]) + 1.0/60.0*(Sour[1][j]-Sour[Ny-6][j]);
		Dest[Ny-2][j] = 0.75*(Sour[0][j]-Sour[Ny-3][j]) - 0.15*(Sour[1][j]-Sour[Ny-4][j]) + 1.0/60.0*(Sour[2][j]-Sour[Ny-5][j]);
		Dest[Ny-1][j] = Dest[0][j];
		//
		for (int i=0; i<Ny; ++i)		Dest[i][j] *= 1.0/dY;
	}
}



void GD2DParXX(const matrix<double>& Sour,matrix<double>& Dest,double dX)
{
#pragma omp parallel for
	for (int i=0; i<Sour.m_Row; ++i) {
		for (int j = 3; j < Sour.m_Column-3; ++j) {
			Dest[i][j] = -49.0/18.0*Sour[i][j] + 1.5*(Sour[i][j+1]+Sour[i][j-1]) - 0.15*(Sour[i][j+2]+Sour[i][j-2]) + 1.0/90.0*(Sour[i][j+3]+Sour[i][j-3]);
		}
		Dest[i][2] = -49.0/18.0*Sour[i][2] + 1.5*(Sour[i][3]+Sour[i][1]) - 0.15*(Sour[i][4]+Sour[i][0]) + 1.0/90.0*(Sour[i][5]+Sour[i][Sour.m_Column-2]);
		Dest[i][1] = -49.0/18.0*Sour[i][1] + 1.5*(Sour[i][2]+Sour[i][0]) - 0.15*(Sour[i][3]+Sour[i][Sour.m_Column-2]) + 1.0/90.0*(Sour[i][4]+Sour[i][Sour.m_Column-3]);
		Dest[i][0] = -49.0/18.0*Sour[i][0] + 1.5*(Sour[i][1]+Sour[i][Sour.m_Column-2]) - 0.15*(Sour[i][2]+Sour[i][Sour.m_Column-3]) + 1.0/90.0*(Sour[i][3]+Sour[i][Sour.m_Column-4]);
		Dest[i][Sour.m_Column-3] = -49.0/18.0*Sour[i][Sour.m_Column-3] + 1.5*(Sour[i][Sour.m_Column-2]+Sour[i][Sour.m_Column-4]) - 0.15*(Sour[i][0]+Sour[i][Sour.m_Column-5]) + 1.0/90.0*(Sour[i][1]+Sour[i][Sour.m_Column-6]);
		Dest[i][Sour.m_Column-2] = -49.0/18.0*Sour[i][Sour.m_Column-2] + 1.5*(Sour[i][0]+Sour[i][Sour.m_Column-3]) - 0.15*(Sour[i][1]+Sour[i][Sour.m_Column-4]) + 1.0/90.0*(Sour[i][2]+Sour[i][Sour.m_Column-5]);
		Dest[i][Sour.m_Column-1] = Dest[i][0];
		//
		for (int j=0; j<Sour.m_Column; ++j)		Dest[i][j] *= 1.0/(dX*dX);
	}
}


void GD2DParYY(const matrix<double>& Sour,matrix<double>& Dest,double dY)
{
	int Nx = Sour.m_Column;
	int Ny = Sour.m_Row;
#pragma omp parallel for
	for (int j=0; j<Nx; ++j) {
		for (int i=3; i<Ny-3; ++i) {
			Dest[i][j] = -49.0/18.0*Sour[i][j] + 1.5*(Sour[i+1][j]+Sour[i-1][j]) - 0.15*(Sour[i+2][j]+Sour[i-2][j]) + 1.0/90.0*(Sour[i+3][j]+Sour[i-3][j]);
		}
		Dest[2][j] = -49.0/18.0*Sour[2][j] + 1.5*(Sour[3][j]+Sour[1][j]) - 0.15*(Sour[4][j]+Sour[0][j]) + 1.0/90.0*(Sour[5][j]+Sour[Ny-2][j]);
		Dest[1][j] = -49.0/18.0*Sour[1][j] + 1.5*(Sour[2][j]+Sour[0][j]) - 0.15*(Sour[3][j]+Sour[Ny-2][j]) + 1.0/90.0*(Sour[4][j]+Sour[Ny-3][j]);
		Dest[0][j] = -49.0/18.0*Sour[0][j] + 1.5*(Sour[1][j]+Sour[Ny-2][j]) - 0.15*(Sour[2][j]+Sour[Ny-3][j]) + 1.0/90.0*(Sour[3][j]+Sour[Ny-4][j]);
		Dest[Ny-3][j] = -49.0/18.0*Sour[Ny-3][j] + 1.5*(Sour[Ny-2][j]+Sour[Ny-4][j]) - 0.15*(Sour[0][j]+Sour[Ny-5][j]) + 1.0/90.0*(Sour[1][j]+Sour[Ny-6][j]);
		Dest[Ny-2][j] = -49.0/18.0*Sour[Ny-2][j] + 1.5*(Sour[0][j]+Sour[Ny-3][j]) - 0.15*(Sour[1][j]+Sour[Ny-4][j]) + 1.0/90.0*(Sour[2][j]+Sour[Ny-5][j]);
		Dest[Ny-1][j] = Dest[0][j];
		//
		for (int i=0; i<Ny; ++i)	Dest[i][j] *= 1.0/(dY*dY);
	}
}


void GD2DLaplace(const matrix<double>& Sour,matrix<double>& Dest,double dX,double dY)
{
	int Nx = Sour.m_Column;
	int Ny = Sour.m_Row;
#pragma omp parallel for
	for (int i=0; i<Sour.m_Row; ++i) {
		for (int j = 3; j < Sour.m_Column-3; ++j) {
			Dest[i][j] = -49.0/18.0*Sour[i][j] + 1.5*(Sour[i][j+1]+Sour[i][j-1]) - 0.15*(Sour[i][j+2]+Sour[i][j-2]) + 1.0/90.0*(Sour[i][j+3]+Sour[i][j-3]);
		}
		Dest[i][2] = -49.0/18.0*Sour[i][2] + 1.5*(Sour[i][3]+Sour[i][1]) - 0.15*(Sour[i][4]+Sour[i][0]) + 1.0/90.0*(Sour[i][5]+Sour[i][Sour.m_Column-2]);
		Dest[i][1] = -49.0/18.0*Sour[i][1] + 1.5*(Sour[i][2]+Sour[i][0]) - 0.15*(Sour[i][3]+Sour[i][Sour.m_Column-2]) + 1.0/90.0*(Sour[i][4]+Sour[i][Sour.m_Column-3]);
		Dest[i][0] = -49.0/18.0*Sour[i][0] + 1.5*(Sour[i][1]+Sour[i][Sour.m_Column-2]) - 0.15*(Sour[i][2]+Sour[i][Sour.m_Column-3]) + 1.0/90.0*(Sour[i][3]+Sour[i][Sour.m_Column-4]);
		Dest[i][Sour.m_Column-3] = -49.0/18.0*Sour[i][Sour.m_Column-3] + 1.5*(Sour[i][Sour.m_Column-2]+Sour[i][Sour.m_Column-4]) - 0.15*(Sour[i][0]+Sour[i][Sour.m_Column-5]) + 1.0/90.0*(Sour[i][1]+Sour[i][Sour.m_Column-6]);
		Dest[i][Sour.m_Column-2] = -49.0/18.0*Sour[i][Sour.m_Column-2] + 1.5*(Sour[i][0]+Sour[i][Sour.m_Column-3]) - 0.15*(Sour[i][1]+Sour[i][Sour.m_Column-4]) + 1.0/90.0*(Sour[i][2]+Sour[i][Sour.m_Column-5]);
		Dest[i][Sour.m_Column-1] = Dest[i][0];
		//
		for (int j=0; j<Sour.m_Column; ++j)		Dest[i][j] *= 1.0/(dX*dX);
	}

#pragma omp parallel for
	for (int j=0; j<Nx; ++j) {
		for (int i=3; i<Ny-3; ++i) {
			Dest[i][j] += 1.0/(dY*dY)*(-49.0/18.0*Sour[i][j] + 1.5*(Sour[i+1][j]+Sour[i-1][j]) - 0.15*(Sour[i+2][j]+Sour[i-2][j]) + 1.0/90.0*(Sour[i+3][j]+Sour[i-3][j]));
		}
		Dest[2][j] += 1.0/(dY*dY)*(-49.0/18.0*Sour[2][j] + 1.5*(Sour[3][j]+Sour[1][j]) - 0.15*(Sour[4][j]+Sour[0][j]) + 1.0/90.0*(Sour[5][j]+Sour[Ny-2][j]));
		Dest[1][j] += 1.0/(dY*dY)*(-49.0/18.0*Sour[1][j] + 1.5*(Sour[2][j]+Sour[0][j]) - 0.15*(Sour[3][j]+Sour[Ny-2][j]) + 1.0/90.0*(Sour[4][j]+Sour[Ny-3][j]));
		Dest[0][j] += 1.0/(dY*dY)*(-49.0/18.0*Sour[0][j] + 1.5*(Sour[1][j]+Sour[Ny-2][j]) - 0.15*(Sour[2][j]+Sour[Ny-3][j]) + 1.0/90.0*(Sour[3][j]+Sour[Ny-4][j]));
		Dest[Ny-3][j] += 1.0/(dY*dY)*(-49.0/18.0*Sour[Ny-3][j] + 1.5*(Sour[Ny-2][j]+Sour[Ny-4][j]) - 0.15*(Sour[0][j]+Sour[Ny-5][j]) + 1.0/90.0*(Sour[1][j]+Sour[Ny-6][j]));
		Dest[Ny-2][j] += 1.0/(dY*dY)*(-49.0/18.0*Sour[Ny-2][j] + 1.5*(Sour[0][j]+Sour[Ny-3][j]) - 0.15*(Sour[1][j]+Sour[Ny-4][j]) + 1.0/90.0*(Sour[2][j]+Sour[Ny-5][j]));
		Dest[Ny-1][j] += Dest[0][j];
		//
		//for (int i=0; i<Ny; ++i)	Dest[i][j] *= 1.0/(dY*dY);
	}
}
