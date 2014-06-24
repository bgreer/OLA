#include "header.h"

// this file contains the function evaluation
// so important it deserves its own file

// chi = sum_space[ (model - target)^2. ]
// model is the linear combination of kernels, coeffients given by parameters
// target is specified and constant
// sum over space is for every point in the solution grid

/*
float feval(int numparams, float *params, measurement *list, float *target, kernel_desc *allkers)
{
	float model;
	float ret = 0.0;
	int ii, iz;

	for (iz=0; iz<allkers->nz; iz++)
	{
		model = 0.0;
		for (ii=0; ii<numparams; ii++)
		{
			// add kernel[ii].sens[iz] * param[ii] to function, subtract from target[iz]
			model += params[ii] * list[ii].ker->sens[0][iz];
		}
		ret += pow( model - target[iz] ,2.0)/(target[iz]+0.1);
	}
	for (ii=0; ii<numparams; ii++)
	{
		ret += 1e-6 * pow(params[ii] * list[ii].ey, 2.0);
	}
	return ret;
}
*/

// reconstruct the velocity and associated error based on optimized parameters
void velocity(int numparams, float *params, measurement *list, float *vx, float *ex, float *vy, float *ey)
{
	int ii, ij;
	*vx = 0.0;
	*ex = 0.0;
	*vy = 0.0;
	*ey = 0.0;

	for (ii=0; ii<numparams; ii++)
	{
		*vx += params[ii] * list[ii].vx;
		*ex += pow(params[ii] * list[ii].ex, 2.0);
		*vy += params[ii] * list[ii].vy;
		*ey += pow(params[ii] * list[ii].ey, 2.0);
	}
	*ex = sqrt(*ex);
	*ey = sqrt(*ey);
}
