#include "header.h"

void precompute_interp (interp_mask *mask)
{

}


/* bicubic interpolation */
void interp (float x, float y, float ***data, int dim, int nz, float *res)
{
	int ii, ij, iz;
	int ix[4], iy[4];
	float r[4][nz], weight, temp;

	if (x < 0.0 || x >= dim || y < 0.0 || y >= dim)
	{
		memset(res, 0x00, nz * sizeof(float));
		return;
	}

	// find nearest pixels
	ix[1] = imax(floor(x), 0);
	iy[1] = imax(floor(y), 0);
	ix[2] = imin(ix[1]+1, dim-1);
	iy[2] = imin(iy[1]+1, dim-1);
	ix[0] = imax(ix[1]-1,0);
	iy[0] = imax(iy[1]-1,0);
	ix[3] = imin(ix[1]+2,dim-1);
	iy[3] = imin(iy[1]+2,dim-1);


	// interpolate at current y
	for (ii=0; ii<4; ii++) // x index
	{
		memset(r[ii], 0x00, nz*sizeof(float));
		weight = 0.0;
		for (ij=0; ij<4; ij++) // y index
		{
			temp = interp_kernel(iy[ij] - y);
			weight += temp;
			for (iz=0; iz<nz; iz++)
				r[ii][iz] += temp * data[ix[ii]][iy[ij]][iz];
		}
		for (iz=0; iz<nz; iz++)
			r[ii][iz] /= weight; // normalize
	}

	memset(res, 0x00, nz*sizeof(float));
	// interpolate at current x
	for (ii=0; ii<4; ii++)
	{
		temp = interp_kernel(ix[ii]-x);
		weight += temp;
		for (iz=0; iz<nz; iz++)
			res[iz] += temp * r[ii][iz];
	}

	for (iz=0; iz<nz; iz++)
		res[iz] *= 2.0/weight;
}

/* for bicubic interpolation */
float interp_kernel (float x0)
{
	float x;
	x = fabs(x0);
	if (x <= 1.0)
		return 1.0 + x*x*(-2.5 + 1.5*x);
	else if (x < 2.0)
		return 2.0 + x*(-4.0 + x*(2.5 - 0.5*x));
	else
		return 0.0;
}

int imin (int a, int b)
{
	return (a > b) ? b : a;
}

int imax (int a, int b)
{
	return (a > b) ? a : b;
}
