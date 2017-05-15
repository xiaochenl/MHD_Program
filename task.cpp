#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif
#include "Task.h"


void Test_Advection_WENO5(void)
{
	const double pi = 3.14159265358979323846;
	double xInf = -1.0;
	double xSup = 1.0;
	double Speed = 2.0;
	double visc = 0.000;
	double timeInterval = 0.0;
	double timeRemain = 0.0;

	int Ngrid = 641;

	Grid1D GD(xInf,xSup,Ngrid);
	CSV1D A(GD);
	A.visc = visc;


	double *Func0 = new double [Ngrid];

	double L = xSup-xInf;
	double k = 2.0*pi/L    *2;
	double T0 = 0.1*50.0*L/Speed;
	double cfl = 0.48;
	double ErrInf = 0.0;

	double RT_start = 0.0;
	double RT_end = 0.0;

	timeInterval = cfl*GD.dX/Speed;
	timeInterval = min(timeInterval,pow(GD.dX,5.0/3.0));
	int Nstep = floor(T0/timeInterval)+1;
	timeRemain = T0/Nstep;

#pragma omp parallel for
	for (int _I = 0; _I < Ngrid; ++_I) {
		A.Data[_I] = sin(k*GD.XG[_I]);
		//B.Func(1,_I+1) = A.Func[_I];
		Func0[_I] = sin(k*(GD.XG[_I]-Speed*T0))*exp(-visc*k*k*T0);
	}
	A.Data[Ngrid-1] = A.Data[0];
	//B.Func(1,Ngrid) = A.Func[0];
	Func0[Ngrid-1] = Func0[0];


	cout << "Time steps is " << Nstep << endl;
	cout << "Time interval is " << timeRemain << endl;
	//bool reach = false;
	RT_start = omp_get_wtime();
	for (int _T = 1; _T <= Nstep; ++_T) {
		A.TimeAdvanceTEST(timeRemain,timeInterval,cfl,Speed);
		//cout << reach << endl;
	}
	RT_end = omp_get_wtime();
	for (int _I = 0; _I < Ngrid; ++_I) {
		if ( fabs(A.Data[_I]-Func0[_I])>ErrInf )	ErrInf =  fabs(A.Data[_I]-Func0[_I]);
	}



	cout << "Max ERROR is " << ErrInf << endl;
	cout << "Run time is " << RT_end-RT_start << endl;
	///////////////////////////////////////////////////////////////////
//	RT_start = omp_get_wtime();
//	for (int _T = 1; _T <= Nstep; ++_T) {
//		B.TimeAdvance(timeRemain,timeInterval,cfl);
//		//cout << "reach" << endl;
//	}
//	RT_end = omp_get_wtime();
//	for (int _I = 0; _I < Ngrid; ++_I) {
//		//cout << "I = " << _I << endl;
//		if ( fabs(B.Func(1,_I+1)-Func0[_I])>ErrInf )	ErrInf =  fabs(B.Func(1,_I+1)-Func0[_I]);
//		//cout << "Error = " << ErrInf << endl;
//	}
//
//
//
//	cout << "Max ERROR is " << ErrInf << endl;
//	cout << "Run time is " << RT_end-RT_start << endl;


	delete [] Func0;
	Func0 = 0;
}



void TestAdvection2D(void)
{
	const double pi = 3.1415926535897932384626;
	double xInf = -0.5,xSup = 0.5,yInf = -0.5,ySup = 0.5;
	int NGX = 201,NGY=201;
	double Vx = 1.0,Vy = 1.0;

	Grid2D GD(xInf,xSup,yInf,ySup,NGX,NGY);
	CSV2D A(GD);
	matrix<double> Func(NGY,NGX);

	double T0 = 10.0;
	double cfl = 0.8;
	double ErrInf = 0.0;

	double RT_start = 0.0;
	double RT_end = 0.0;
	double timeInterval;
	//double timeRemain;

	double x,y;
#pragma omp parallel for private(x,y)
	for (int i=0; i<NGX; ++i) {
		for (int j=0; j<NGY; ++j) {
			x = GD.XG[i];
			y = GD.YG[j];
			A.Data[j][i] = sin(2.0*pi*x+2.0*pi*y);
			Func[j][i] = sin(2.0*pi*(x-Vx*T0)+2.0*pi*(y-Vy*T0));
		}
	}

	bool reach = false;
	double currentTime = 0.0;
	RT_start = omp_get_wtime();
	while (!reach) {
		reach = A.TimeAdvanceTEST(T0-currentTime,timeInterval,cfl,Vx,Vy);
		currentTime += timeInterval;
		cout << "Current time is " << currentTime << endl;
	}
	RT_end = omp_get_wtime();

	for (int i=0; i<NGX; ++i) {
		for (int j=0; j<NGY; ++j) {
			if ( fabs(A.Data[j][i]-Func[j][i])>ErrInf )	ErrInf =  fabs(A.Data[j][i]-Func[j][i]);
		}
	}



	cout << "Max ERROR is " << ErrInf << endl;
	cout << "Run time is " << RT_end-RT_start << endl;
}


void GLM_MHD2D_Run(void)
{
	//string Sample = "OrszagTangVortex";
	//string Sample = "CoalescenceInstability";
	//string Sample = "Boom";
	string Sample = "MagneticCloud";

	const double pi = 3.1415926535897932384626;
	double xInf,xSup,yInf,ySup;
	int Nx,Ny,outputFrames;
	double startTime,endTime;

	
	if (Sample=="OrszagTangVortex") {
		xInf = 0.0;
		xSup = 2.0*pi;
		yInf = 0.0;
		ySup = 2.0*pi;

	} else if (Sample=="CoalescenceInstability") {
		xInf = -0.5;
		xSup = 0.5;
		yInf = -0.5;
		ySup = 0.5;

	}
	  else if (Sample == "Boom") {
		xInf = -0.5;
		xSup = 0.5;
		yInf = -0.5;
		ySup = 0.5;

	}
	  else if (Sample == "MagneticCloud") {
		  xInf = -1.0;
		  xSup = 1.0;
		  yInf = -1.0;
		  ySup = 1.0;

	  }


	Nx = 151;
	Ny = 151;

	matrix<double> XXG(1,Nx),YYG(1,Ny);

	Grid2D GD(xInf,xSup,yInf,ySup,Nx,Ny);
	GLM_MHD2D A(GD);


	if (Sample=="OrszagTangVortex") {
		OrszagTangVortex(A);
	} else if (Sample=="CoalescenceInstability") {
		CoalescenceInstability(A);
	}
	  else if (Sample == "Boom") {
		Boom(A);
	}
	  else if (Sample == "MagneticCloud") {
		  MagneticCloud(A);
	  }
	A.Prim2CSV();

	/////////////////////////////////////////////
#pragma omp parallel
	{
#pragma omp for
		for (int i=0; i<Nx; ++i) {
			XXG[0][i] = GD.XG[i];
		}
#pragma omp for
		for (int i=0; i<Ny; ++i) {
			YYG[0][i] = GD.YG[i];
		}
	}
	XXG.WriteToFile("E:/output/X");
	YYG.WriteToFile("E:/output/Y");
	//////////////////////////////////////////////

	bool reach;
	outputFrames = 1;
	startTime = 0.0; endTime = 0.1;
	double xInterval = (xSup-xInf)/(Nx-1);
	double yInterval = (ySup-yInf)/(Ny-1);
	double currentTime = startTime;
	matrix<double> timePoints(1,outputFrames),EK(1,outputFrames),EB(1,outputFrames),EP(1,outputFrames),ET(1,outputFrames);
	double dS = xInterval*yInterval;


//	for (int i=0;i<outputFrames;++i)
//	{
//		timePoints[0][i] = startTime + (i+1)*(endTime-startTime)/outputFrames;
//
//	}

	
	timePoints[0][0] = endTime;
	//timePoints[0][1] = 11.5;
	//timePoints[0][2] = 30.0;

	//timePoints.WriteToFile("E:/output/T");

	int targetFrame = 0;
	double timeInterval;
	double cfl = 0.6;


	//	Start time evolution

	while (targetFrame < outputFrames ) {
		reach = A.TimeAdvance(timePoints[0][targetFrame]-currentTime,timeInterval,cfl);
		if (!reach)
		{
		currentTime += timeInterval;
		cout << "Current time is " << currentTime << ", time interval is " << timeInterval << endl;
		}
		else
		{
			EK[0][targetFrame] = 0.0;
			EB[0][targetFrame] = 0.0;
			ET[0][targetFrame] = 0.0;
			EP[0][targetFrame] = 0.0;

			for (int i = 0;i<Ny-1;++i)
			{
				for (int j=1; j<Nx-1;++j)
				{
				EK[0][targetFrame] += 0.5*dS*(A.Vx[i][j]*A.Vx[i][j]+A.Vy[i][j]*A.Vy[i][j]+A.Vz[i][j]*A.Vz[i][j]);
				EB[0][targetFrame] += 0.5*dS*(A.Bx[i][j]*A.Bx[i][j]+A.By[i][j]*A.By[i][j]+A.Bz[i][j]*A.Bz[i][j]);
				ET[0][targetFrame] += dS*A.Ener[i][j];
				}
			}

			EP[0][targetFrame] = ET[0][targetFrame]-EB[0][targetFrame]-EK[0][targetFrame];

			//A.Vx.WriteToFile("Output//Vx",targetFrame);
			A.Mass.Data.WriteToFile("E:/output/Mass");
			A.Pres.WriteToFile("E:/output//Pre");
			A.Psi.Data.WriteToFile("E:/output/Psi");
			//A.Mass.Data.WriteToFile("E:/output//Mass", targetFrame);
			//A.Pres.WriteToFile("E:/output//Pres", targetFrame);
			//A.Psi.Data.WriteToFile("E:/output//Psi", targetFrame);

			//currentTime = timePoints[0][targetFrame];
			//cout  << "Current time is " << currentTime << ", Frame = " << targetFrame << endl;

			targetFrame++ ;
		}
	}
	EK.WriteToFile("E:/output/EK");
	EB.WriteToFile("E:/output/EB");
	ET.WriteToFile("E:/output/ET");
	EP.WriteToFile("E:/output/EP");

}
