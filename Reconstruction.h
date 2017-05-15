/*
 * Reconstruction.h
 *
 *  Created on: Mar 12, 2016
 *      Author: majun
 */

#ifndef RECONSTRUCTION_H_
#define RECONSTRUCTION_H_


class TMP_WENO
{
public:
	TMP_WENO(void)  {eps = 1.0e-40;	_W = 0.0;}
public:
	double V[3],alpha[3],beta[3],W[3],_W;
	double eps;
};

typedef double (*wenoKIND)(double B2,double B1,double _V,double F1,double F2,TMP_WENO& tmp);


double WENO3(double _B, double _V, double _F, TMP_WENO& tmp);

double WENO5(double B2, double B1, double _V, double F1, double F2, TMP_WENO& tmp);

double WENO5M(double B2, double B1, double _V, double F1, double F2, TMP_WENO& tmp);
#endif /* RECONSTRUCTION_H_ */
