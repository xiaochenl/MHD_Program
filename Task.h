/*
 * Task.h
 *
 *  Created on: Mar 14, 2016
 *      Author: majun
 */

#ifndef TASK_H_
#define TASK_H_

#include <iostream>
#include <cmath>
#include <omp.h>
#include "Matrix.h"
#include "scheme1D.h"
#include "scheme2D.h"
#include "SetInitial.h"
using namespace std;

void Test_Advection_WENO5(void);

void TestAdvection2D(void);

void GLM_MHD2D_Run(void);

#endif /* TASK_H_ */
