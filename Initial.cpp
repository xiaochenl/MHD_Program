/*
 * SetInitial.cpp
 *
 *  Created on: Mar 16, 2016
 *      Author: majun
 */

#include "SetInitial.h"


void OrszagTangVortex(GLM_MHD2D& A)
{
        double gamma = 5.0/3.0;
        double x,y;
#pragma omp parallel for private(x,y)
        for (int i=0; i<(*A.GD).NGY; ++i) {
                for (int j=0; j<(*A.GD).NGX; ++j) {
                        x = (*A.GD).XG[j];
                        y = (*A.GD).YG[i];

                        A.Mass0[i][j] = gamma*gamma*0.0;
                        A.Vx0[i][j] = 0.0;
                        A.Vy0[i][j] = 0.0;
                        A.Vz0[i][j] = 0.0;
                        A.Pres0[i][j] = 0.0*gamma;
                        A.Bx0[i][j] = 0.0;
                        A.By0[i][j] = 0.0;
                        A.Bz0[i][j] = 0.0;

                        A.Mass.Data[i][j] = gamma*gamma*1.0;
                        A.Vx[i][j] = -sin(y);
                        A.Vy[i][j] = sin(x);
                        A.Vz[i][j] = 0.0;
                        A.Pres[i][j] = gamma*1.0;
                        A.Bx.Data[i][j] = -sin(y);
                        A.By.Data[i][j] = sin(2.0*x);
                        A.Bz.Data[i][j] = 0.0;
                        A.Psi.Data[i][j] = 0.0;

//			A.Vx[i][j] = -sin(y);
//			A.Vy[i][j] = sin(x);
//			A.Bx.Data[i][j] = -sin(2.0*y);
//			A.By.Data[i][j] = sin(x);
                }
        }
}

void Boom(GLM_MHD2D& A)
{
        double gamma = 5.0 / 3.0;
        double x, y;
#pragma omp parallel for private(x,y)
        for (int i = 0; i < (*A.GD).NGY; ++i) {
                for (int j = 0; j < (*A.GD).NGX; ++j) {
                        x = (*A.GD).XG[j];
                        y = (*A.GD).YG[i];
                        if (sqrt(x*x) + sqrt(y*y) < 0.1)
                        {
                                A.Mass0[i][j] = gamma*gamma*0.0;
                                A.Vx0[i][j] = 0.0;
                                A.Vy0[i][j] = 0.0;
                                A.Vz0[i][j] = 0.0;
                                A.Pres0[i][j] = 0.0*gamma;
                                A.Bx0[i][j] = 0.0;
                                A.By0[i][j] = 0.0;
                                A.Bz0[i][j] = 0.0;

                                A.Mass.Data[i][j] = 1.0;
                                A.Vx[i][j] = 0;
                                A.Vy[i][j] = 0;
                                A.Vz[i][j] = 0.0;
                                A.Pres[i][j] = 10.0;
                                A.Bx.Data[i][j] = 1 / sqrt(2);
                                A.By.Data[i][j] = 1 / sqrt(2);
                                A.Bz.Data[i][j] = 0.0;
                                A.Psi.Data[i][j] = 0.0;
                        }
                        else
                        {
                                A.Mass0[i][j] = gamma*gamma*0.0;
                                A.Vx0[i][j] = 0.0;
                                A.Vy0[i][j] = 0.0;
                                A.Vz0[i][j] = 0.0;
                                A.Pres0[i][j] = 0.0*gamma;
                                A.Bx0[i][j] = 0.0;
                                A.By0[i][j] = 0.0;
                                A.Bz0[i][j] = 0.0;

                                A.Mass.Data[i][j] = 1.0;
                                A.Vx[i][j] = 0;
                                A.Vy[i][j] = 0;
                                A.Vz[i][j] = 0.0;
                                A.Pres[i][j] = 0.1;
                                A.Bx.Data[i][j] = 1 / sqrt(2);
                                A.By.Data[i][j] = 1 / sqrt(2);
                                A.Bz.Data[i][j] = 0.0;
                                A.Psi.Data[i][j] = 0.0;
                        }
                }
        }
}

void MagneticCloud(GLM_MHD2D& A)
{
        double gamma = 5.0 / 3.0;
        double x, y;
#pragma omp parallel for private(x,y)
        for (int i = 0; i < (*A.GD).NGY; ++i) {
                for (int j = 0; j < (*A.GD).NGX; ++j) {
                        x = (*A.GD).XG[j];
                        y = (*A.GD).YG[i];
                        A.Mass0[i][j] = gamma*gamma*0.0;
                        A.Vx0[i][j] = 0.0;
                        A.Vy0[i][j] = 0.0;
                        A.Vz0[i][j] = 0.0;
                        A.Pres0[i][j] = 0.0*gamma;
                        A.Bx0[i][j] = 0.0;
                        A.By0[i][j] = 0.0;
                        A.Bz0[i][j] = 0.0;

                        if (x >=0.0 && y >= 0.0)
                        {
                                A.Mass.Data[i][j] = 0.9308;
                                A.Vx[i][j] = 1.5639;
                                A.Vy[i][j] = -0.4977;
                                A.Vz[i][j] = 0.0618;
                                A.Pres[i][j] = 5.0838;
                                A.Bx.Data[i][j] = 0.3501;
                                A.By.Data[i][j] = 0.9830;
                                A.Bz.Data[i][j] = 0.3050;
                                A.Psi.Data[i][j] = 0.0;
                        }
                        //%		else if ((x - 0.25)*(x - 0.25) + (y - 0.5)*(y - 0.5)<0.0225){
                        //			A.Mass0[i][j] = gamma*gamma*0.0;
                        //			A.Vx0[i][j] = 0.0;
                        //			A.Vy0[i][j] = 0.0;
                        //			A.Vz0[i][j] = 0.0;
                        //			A.Pres0[i][j] = 0.0*gamma;
                        //			A.Bx0[i][j] = 0.0;
                        //			A.By0[i][j] = 0.0;
                        //			A.Bz0[i][j] = 0.0;
                        //
                        //			A.Mass.Data[i][j] = 10;
                        //			A.Vx[i][j] = 0;
                        //			A.Vy[i][j] = 0;
                        //			A.Vz[i][j] = 0.0;
                        //			A.Pres[i][j] = 0;
                        //			A.Bx.Data[i][j] = 0;
                        //			A.By.Data[i][j] = 0;
                        //			A.Bz.Data[i][j] = 0;
                        //			A.Psi.Data[i][j] = 0.0;
                        //		}

                        else if (x < 0.0 && y >= 0.0)
                        {
                                A.Mass.Data[i][j] = 1.0304;
                                A.Vx[i][j] = 1.5309;
                                A.Vy[i][j] = -1.0147;
                                A.Vz[i][j] = -0.0986;
                                A.Pres[i][j] = 5.7813;
                                A.Bx.Data[i][j] = 0.3501;
                                A.By.Data[i][j] = 0.5078;
                                A.Bz.Data[i][j] = 0.1576;
                                A.Psi.Data[i][j] = 0.0;
                        }
                        else if (x < 0.0 && y < 0.0)
                        {
                                A.Mass.Data[i][j] = 1.0;
                                A.Vx[i][j] = 1.7500;
                                A.Vy[i][j] = -1.0000;
                                A.Vz[i][j] = 0.0000;
                                A.Pres[i][j] = 6.0000;
                                A.Bx.Data[i][j] = 0.5642;
                                A.By.Data[i][j] = 0.5078;
                                A.Bz.Data[i][j] = 0.2539;
                                A.Psi.Data[i][j] = 0.0;
                        }
                        else if (x >= 0 && y < 0)
                        {
                                A.Mass.Data[i][j] = 1.8887;
                                A.Vx[i][j] = 0.1236;
                                A.Vy[i][j] = -0.9224;
                                A.Vz[i][j] = 0.0388;
                                A.Pres[i][j] = 12.9990;
                                A.Bx.Data[i][j] = 0.5642;
                                A.By.Data[i][j] = 0.9830;
                                A.Bz.Data[i][j] = 0.4915;
                                A.Psi.Data[i][j] = 0.0;
                        }
                }
        }
}


void CoalescenceInstability(GLM_MHD2D& A)
{
        const double pi = 3.1415926535897932384626;
        //double gamma = 5.0/3.0;
        double x,y,A0;
#pragma omp parallel for private(x,y,A0)
        for (int i=0; i<(*A.GD).NGY; ++i) {
                for (int j=0; j<(*A.GD).NGX; ++j) {
                        x = (*A.GD).XG[j];
                        y = (*A.GD).YG[i];

                        A0 = 0.05*( cos(4.0*pi*x)-cos(4.0*pi*y) );

                        A.Mass0[i][j] = 1.0;
                        A.Vx0[i][j] = 0.0;
                        A.Vy0[i][j] = 0.0;
                        A.Vz0[i][j] = 0.0;
                        A.Pres0[i][j] = 1.0+8.0*pi*pi*A0*A0;
                        A.Bx0[i][j] = 0.2*pi*sin(4.0*pi*y);
                        A.By0[i][j] = 0.2*pi*sin(4.0*pi*x);
                        A.Bz0[i][j] = 0.0;

                        A.Mass.Data[i][j] = 0.0;
                        A.Vx[i][j] = 1.0e-4*sin(2.0*pi*x)*cos(2.0*pi*y);
                        A.Vy[i][j] = -1.0e-4*cos(2.0*pi*x)*sin(2.0*pi*y);
                        A.Vz[i][j] = 0.0;
                        A.Pres[i][j] = 0.0;
                        A.Bx.Data[i][j] = 0.0;
                        A.By.Data[i][j] = 0.0;
                        A.Bz.Data[i][j] = 0.0;
                        A.Psi.Data[i][j] = 0.0;

                }
        }
}

