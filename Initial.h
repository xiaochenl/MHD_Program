/*
 * SetInitial.h
 *
 *  Created on: Mar 16, 2016
 *      Author: majun
 */

#ifndef SETINITIAL_H_
#define SETINITIAL_H_

#include "scheme2D.h"
#include <omp.h>
#include <cmath>

#include "Matrix.h"

void OrszagTangVortex(GLM_MHD2D& A);

void CoalescenceInstability(GLM_MHD2D& A);
void Boom(GLM_MHD2D& A);
void MagneticCloud(GLM_MHD2D& A);

#endif /* SETINITIAL_H_ */
