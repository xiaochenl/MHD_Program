/*
 * Reconstruction.cpp
 *
 *  Created on: Mar 12, 2016
 *      Author: majun
 */



#include "Reconstruction.h"


double WENO3(double _B, double _V, double _F, TMP_WENO& tmp)
{
	tmp.eps = 1.0e-6;

	tmp.V[0] = 0.5*_V + 0.5*_F;
	tmp.V[1] = -0.5*_B + 1.5*_F;

	tmp.beta[0] = (_F-_V)*(_F-_V);
	tmp.beta[1] = (_V-_B)*(_V-_B);

	tmp.alpha[0] = 2.0/3.0/((tmp.beta[0]+tmp.eps)*(tmp.beta[0]+tmp.eps));
	tmp.alpha[1] = 1.0/3.0/((tmp.beta[1]+tmp.eps)*(tmp.beta[1]+tmp.eps));

	tmp._W = tmp.alpha[0]+tmp.alpha[1];
	tmp.W[0] = tmp.alpha[0]/tmp._W;
	tmp.W[1] = tmp.alpha[1]/tmp._W;


	return tmp.W[0]*tmp.V[0] + tmp.W[1]*tmp.V[1];
}




double WENO5(double B2, double B1, double _V, double F1, double F2, TMP_WENO& tmp)
{
	tmp.eps = 1.0e-6;

	tmp.V[0] = 1.0/3.0*_V+5.0/6.0*F1-1.0/6.0*F2;
	tmp.V[1] = -1.0/6.0*B1+5.0/6.0*_V+1.0/3.0*F1;
	tmp.V[2] = 1.0/3.0*B2-7.0/6.0*B1+11.0/6.0*_V;

	tmp.beta[0] = 13.0/12.0*(_V-2.0*F1+F2)*(_V-2.0*F1+F2)+0.25*(3.0*_V-4.0*F1+F2)*(3.0*_V-4.0*F1+F2);
	tmp.beta[1] = 13.0/12.0*(B1-2.0*_V+F1)*(B1-2.0*_V+F1)+0.25*(B1-F1)*(B1-F1);
	tmp.beta[2] = 13.0/12.0*(B2-2.0*B1+_V)*(_V-2.0*B1+B2)+0.25*(3.0*_V-4.0*B1+B2)*(3.0*_V-4.0*B1+B2);

	tmp.alpha[0] = 0.3/((tmp.eps+tmp.beta[0])*(tmp.eps+tmp.beta[0]));
	tmp.alpha[1] = 0.6/((tmp.eps+tmp.beta[1])*(tmp.eps+tmp.beta[1]));
	tmp.alpha[2] = 0.1/((tmp.eps+tmp.beta[2])*(tmp.eps+tmp.beta[2]));

	tmp._W = tmp.alpha[0]+tmp.alpha[1]+tmp.alpha[2];
	tmp.W[0] = tmp.alpha[0]/tmp._W;
	tmp.W[1] = tmp.alpha[1]/tmp._W;
	tmp.W[2] = tmp.alpha[2]/tmp._W;

	return tmp.W[0]*tmp.V[0]+tmp.W[1]*tmp.V[1]+tmp.W[2]*tmp.V[2];
}




double WENO5M(double B2, double B1, double _V, double F1, double F2, TMP_WENO& tmp)
{
	//tmp.eps = 1.0e-6;

	tmp.V[0] = 1.0/3.0*_V+5.0/6.0*F1-1.0/6.0*F2;
	tmp.V[1] = -1.0/6.0*B1+5.0/6.0*_V+1.0/3.0*F1;
	tmp.V[2] = 1.0/3.0*B2-7.0/6.0*B1+11.0/6.0*_V;

	tmp.beta[0] = 13.0/12.0*(_V-2.0*F1+F2)*(_V-2.0*F1+F2)+0.25*(3.0*_V-4.0*F1+F2)*(3.0*_V-4.0*F1+F2);
	tmp.beta[1] = 13.0/12.0*(B1-2.0*_V+F1)*(B1-2.0*_V+F1)+0.25*(B1-F1)*(B1-F1);
	tmp.beta[2] = 13.0/12.0*(B2-2.0*B1+_V)*(_V-2.0*B1+B2)+0.25*(3.0*_V-4.0*B1+B2)*(3.0*_V-4.0*B1+B2);

	tmp.alpha[0] = 0.3/((tmp.eps+tmp.beta[0])*(tmp.eps+tmp.beta[0]));
	tmp.alpha[1] = 0.6/((tmp.eps+tmp.beta[1])*(tmp.eps+tmp.beta[1]));
	tmp.alpha[2] = 0.1/((tmp.eps+tmp.beta[2])*(tmp.eps+tmp.beta[2]));

	tmp._W = tmp.alpha[0]+tmp.alpha[1]+tmp.alpha[2];
	tmp.W[0] = tmp.alpha[0]/tmp._W;
	tmp.W[1] = tmp.alpha[1]/tmp._W;
	tmp.W[2] = tmp.alpha[2]/tmp._W;

	//////////////////////////////////////////////////
	tmp.alpha[0] = tmp.W[0]*(0.39-0.9*tmp.W[0]+tmp.W[0]*tmp.W[0])/(0.09+0.4*tmp.W[0]);
	tmp.alpha[1] = tmp.W[1]*(0.96-1.8*tmp.W[1]+tmp.W[1]*tmp.W[1])/(0.36-0.2*tmp.W[1]);
	tmp.alpha[2] = tmp.W[2]*(0.11-0.3*tmp.W[2]+tmp.W[2]*tmp.W[2])/(0.01+0.8*tmp.W[2]);

	tmp._W = tmp.alpha[0]+tmp.alpha[1]+tmp.alpha[2];
	tmp.W[0] = tmp.alpha[0]/tmp._W;
	tmp.W[1] = tmp.alpha[1]/tmp._W;
	tmp.W[2] = tmp.alpha[2]/tmp._W;

	return tmp.W[0]*tmp.V[0]+tmp.W[1]*tmp.V[1]+tmp.W[2]*tmp.V[2];
}
