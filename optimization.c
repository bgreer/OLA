#include "header.h"
#include "/home8/begr7169/SOFTWARE/lapack-3.4.0/lapacke/include/lapacke.h"

/* problem is too big for matrix inversion, use gradient descent on blocks of parameters 
 * possible blocking schemes:
 *  - random blocks
 *  - proximity to target
 *  - dimensions
 *
 *  need to give it the number of parameters to use in the block
 *  an array of which measurements to include in K matrix (length=numparams)
 *  an array of all parameters, we will change the appropriate ones based on the previous array
 *  how much to regularize
 *  what the target function looks like
 *
 *  numparams = number of parameters to invert for
 *  ind[] = array of indices, length numparams, giving which params from main array to use
 *  params[] = full array of parameters, we will alter some of them (params[ind[*]])
 *  mlist[] = full array of measurements, they contain pointers to kernels
 *  target[][][] = target function
 *  allkers = struct containing useful info related to all kernels
 *  lambda = regularization value
*/
void block_opt (int numparams, int *ind, float *params, measurement *mlist, target *t, avgker *avg, grid *g, kernel_desc *allkers, int numparams_all, measurement *****mgrid)
{
	int ii, ij, ix, iy, iz, ki, kj, index, ik;
	int x0, y0, x1, y1, x2, y2, x3, y3;
	float dist, x, y, val, *res;
	lapack_int n, lda, info, worksize;
	double *ktk, *rscale, *cscale, *work;
	double *rhs, *ans;
	int *piv;
	double row, col, amax, lambda2;

	lambda2 = pow(10.0, -t->lambda);

	ktk = (double*) malloc(numparams * numparams * sizeof(double));
	piv = malloc(numparams*numparams*sizeof(int));
	rscale = malloc(numparams*sizeof(double));
	cscale = malloc(numparams*sizeof(double));

	// grid of [x][y][z] is [x*ny*nz + y*nz + z]

	// load K^T K + R
//	printf("Constructing Matrix..\n");
	for (ii=0; ii<numparams; ii++)
	{
		ki = ind[ii]; // global indices of measurements in question
		
		for (ij=0; ij<numparams; ij++)
		{
			kj = ind[ij];
			
			index = ii*numparams+ij; // index of ktk array
			// init to zero
			ktk[index] = 0.0;

			// if they are at the same position, use pre-computed value
			if (mlist[ki].locx == mlist[kj].locx && mlist[ki].locy == mlist[kj].locy)
			{
				ktk[index] += allkers->overlap
					[mlist[ki].n*allkers->maxk + mlist[ki].k]
					[mlist[kj].n*allkers->maxk + mlist[kj].k];
			} else {
				// find pixel bounds of first kernel
				x0 = (mlist[ki].flocx - 0.5*mlist[ki].size - g->x[0])/g->dx;
				y0 = (mlist[ki].flocy - 0.5*mlist[ki].size - g->y[0])/g->dx;
				x1 = x0 + mlist[ki].size / g->dx;
				y1 = y0 + mlist[ki].size / g->dx;

				x2 = (mlist[kj].flocx - 0.5*mlist[kj].size - g->x[0])/g->dx;
				y2 = (mlist[kj].flocy - 0.5*mlist[kj].size - g->y[0])/g->dx;
				x3 = x2 + mlist[kj].size / g->dx;
				y3 = y2 + mlist[kj].size / g->dx;

				// determine overlap pixel range
				x0 = (x0>x2)?x0:x2; // max
				y0 = (y0>y2)?y0:y2; // max
				x1 = (x1>x3)?x3:x1; // min
				y1 = (y1>y3)?y3:y1; // min
				// check bounds
				if (x0 < 0) x0 = 0;
				if (y0 < 0) y0 = 0;
				if (x1 >= g->nx) x1 = g->nx-1;
				if (y1 >= g->ny) y1 = g->ny-1;


				// determine if these kernels overlap in space at all
				if (x1-x0 > 0 || y1-y0 > 0)
				{
					// sum kernel[ki] * kernel[kj] over space
					for (ix=x0; ix<=x1; ix++)
					{
						for (iy=y0; iy<=y1; iy++)
						{
							dist = sqrt(pow(mlist[ki].flocx - g->x[ix],2.) + pow(mlist[ki].flocy - g->y[iy],2.));
							if (dist <= 0.5*mlist[ki].size)
							{
							dist = sqrt(pow(mlist[kj].flocx - g->x[ix],2.) + pow(mlist[kj].flocy - g->y[iy],2.));
							if (dist <= 0.5*mlist[kj].size)
							{
							for (iz=0; iz<g->nz; iz++)
							{
								printf("ERROR\n");
								exit(-1);
								ktk[index] += mlist[ki].ker->sens[8][8][iz] * mlist[kj].ker->sens[8][8][iz] * allkers->dz[iz] * g->dxdy;
							}
							}
							}
						}
					}
				}
			
			}
			// regularization
			ktk[index] += lambda2 * covar(mlist+ki, mlist+kj, 1, 1);
		}
	}

	// pre-conditioning
//	printf("Pre-conditioning..\n");
	n = numparams;
	lda = numparams;
	dgeequ_(&n, &n, ktk, &lda, rscale, cscale, &row, &col, &amax, &info);
	free(rscale);
	free(cscale);

	// lu decomposition
//	printf("LU Decomposition..\n");
	n = numparams;
	lda = numparams;
	dgetrf_(&n, &n, ktk, &lda, piv, &info);
	if (info) printf("\tLU decomp status: %d\n", info);
	if (info) return;

	// inversion
//	printf("Inversion..\n");
	worksize = 32*n;
	work = malloc(worksize*sizeof(double));
	dgetri_(&n, ktk, &lda, piv, work, &worksize, &info);
	if (info) printf("\tInversion status: %d\n", info);
	if (info) return;

	// compute K^T (T-A)
//	printf("Computing RHS..\n");
	rhs = (double*) malloc(numparams * sizeof(double));
	memset(rhs, 0x00, numparams * sizeof(double));
	res = (float*) malloc(allkers->nz * sizeof(float));
	for (ii=0; ii<numparams; ii++)
	{
		index = ind[ii];
		x0 = (mlist[index].flocx - 0.5*mlist[index].size - g->x[0])/g->dx;
		y0 = (mlist[index].flocy - 0.5*mlist[index].size - g->y[0])/g->dx;
		x1 = x0 + mlist[index].size / g->dx;
		y1 = y0 + mlist[index].size / g->dx;

		// check bounds
		if (x0 < 0) x0 = 0;
		if (y0 < 0) y0 = 0;
		if (x1 >= g->nx) x1 = g->nx-1;
		if (y1 >= g->ny) y1 = g->ny-1;

		for (ix=x0; ix<=x1; ix++)
		{
			x = ((allkers->nx-1.)/2.) + (g->x[ix]-mlist[ind[ii]].flocx)/allkers->dx;
			for (iy=y0; iy<=y1; iy++)
			{
				y = ((allkers->ny-1.)/2.) + (g->y[iy]-mlist[ind[ii]].flocy)/allkers->dy;
				interp(x, y, mlist[ind[ii]].ker->sens, allkers->nx, allkers->nz, res);
				for (iz=0; iz<g->nz; iz++)
				{
					rhs[ii] += res[iz] * 
							(t->sens[ix][iy][iz] - avg->sens[ix][iy][iz]) * 
							allkers->dz[iz] * g->dxdy;
				}
			}
		}
		// add covariance
		// add covar(this a, all b's that overlap)
		for (ix=-allkers->numoverlap; ix<=allkers->numoverlap; ix++)
		{
			for (iy=-allkers->numoverlap; iy<=allkers->numoverlap; iy++)
			{
				x0 = mlist[ii].locx + ix;
				y0 = mlist[ii].locy + iy;
				if (x0 >= 0 && y0 >= 0 && x0 <= allkers->maxx && y0 <= allkers->maxy)
					for (ki=0; ki<=allkers->maxn; ki++)
						for (kj=0; kj<=allkers->maxk; kj++)
							if (mgrid[x0][y0][ki][kj] != NULL)
								rhs[ii] -= 2.0 * lambda2 * params[mgrid[x0][y0][ki][kj]->myindex] * covar(mlist+ii, mgrid[x0][y0][ki][kj], 1, 0);
						
			}
		}
	}
	free(res);

//	printf("Computing Updated Parameters..\n");
	ans = (double*) malloc(numparams * sizeof(double));
	memset(ans, 0x00, numparams * sizeof(double));
	for (ii=0; ii<numparams; ii++)
		for (ij=0; ij<numparams; ij++)
			ans[ii] += ktk[ii*numparams + ij] * rhs[ij];

	// update params
	for (ii=0; ii<numparams; ii++)
		params[ind[ii]] += ans[ii];

	free(ktk);
	free(ans);
	free(rhs);
	free(work);
	free(piv);
}

// compute the cost function
void compute_cost (int numparams, float *params, measurement *mlist, target *t, avgker *avg, grid *g, kernel_desc *allkers, float *diff, float *reg, float *norm)
{
	int ix, iy, iz, ii, ij;
	float res, err;

	res = 0.0;
	*norm = 0.0;

	for (ix=0; ix<g->nx; ix++)
	{
		for (iy=0; iy<g->ny; iy++)
		{
			for (iz=0; iz<g->nz; iz++)
			{
				res += pow(t->sens[ix][iy][iz] - avg->sens[ix][iy][iz], 2.0) * allkers->dz[iz] * g->dxdy;
				*norm += avg->sens[ix][iy][iz] * allkers->dz[iz] * g->dxdy;
			}
		}
	}
/*
	err = 0.0;
	for (ii=0; ii<numparams; ii++)
	{
		if (params[ii] != 0.0)
		{
			for (ij=0; ij<numparams; ij++)
			{
				if (params[ij] != 0.0)
				{
//					printf("2 %d %f\n", ii, params[ii] * params[ii]);
					err += params[ii] * params[ij] * covar(mlist+ii, mlist+ij, -1, -1);
				}
			}
		}
	}
*/
	*diff = res;
//	*reg = sqrt(err);
	*reg = 0.0;
}




/* solve problem using one big matrix inversion */
void inv_opt(int numparams, float *params, measurement *list, float *target, kernel_desc *allkers, float lambda)
{
	double *ktk, *rscale, *cscale, *work, *ans, *ans2;
	int *piv;
	double row, col, amax, temp;
	lapack_int n, lda, info, worksize;
	int ii, ij, iz;
	
	ktk = (double*) malloc(numparams * numparams * sizeof(double));
	piv = malloc(numparams*numparams*sizeof(int));
	rscale = malloc(numparams*sizeof(double));
	cscale = malloc(numparams*sizeof(double));

	// load ktk
	for (ii=0; ii<numparams; ii++)
	{
		for (ij=0; ij<numparams; ij++)
		{
			ktk[ii*numparams+ij] = 0.0;
			for (iz=0; iz<allkers->nz; iz++)
				ktk[ii*numparams+ij] += list[ii].ker->sens[0][0][iz] * list[ij].ker->sens[0][0][iz];
			// regularization
			if (ii == ij)
				ktk[ii*numparams+ij] += pow(10.,-lambda) * list[ii].ey * list[ii].ey;
//			printf("%d\t%d\t%e\n", ii,ij,ktk[ii*numparams+ij]);
		}
//		printf("\n");
	}

	// pre-conditioning
	n = numparams;
	lda = numparams;
	dgeequ_(&n, &n, ktk, &lda, rscale, cscale, &row, &col, &amax, &info);
	free(rscale);
	free(cscale);

	// lu decomposition
	n = numparams;
	lda = numparams;
	dgetrf_(&n, &n, ktk, &lda, piv, &info);
	if (info) printf("\tLU decomp status: %d\n", info);

	// inversion
	worksize = 32*n;
	work = malloc(worksize*sizeof(double));
	dgetri_(&n, ktk, &lda, piv, work, &worksize, &info);
	if (info) printf("\tInversion status: %d\n", info);

	free(work);
	free(piv);
	ans = (double*) malloc(numparams * sizeof(double));
	ans2 = (double*) malloc(numparams * sizeof(double));

	for (ii=0; ii<numparams; ii++)
	{
		ans[ii] = 0.0;
		for (iz=0; iz<allkers->nz; iz++)
		{
			ans[ii] += list[ii].ker->sens[0][0][iz] * target[iz];
		}
//		printf("%e\n", ans[ii]);
	}
	// ans[ii] now contains K^T * T
	// multiply (ktk^-1) * (K^T * T) to get ans2
	for (ii=0; ii<numparams; ii++)
	{
		ans2[ii] = 0.0;
		for (ij=0; ij<numparams; ij++)
		{
			ans2[ii] += ktk[ii*numparams + ij] * ans[ij];
		}
	}
	free(ktk);

	// copy results out
	for (ii=0; ii<numparams; ii++)
		params[ii] = ans2[ii];
	free(ans);
	free(ans2);
}

/*
void line_opt(int numparams, float *params, measurement *list, float *target, kernel_desc *allkers)
{
	float chi, lastchi, tol;
	float step[numparams];
	int ii, ij, iter, converged;
	int ik;

	tol = 1e-6;
	for (ii=0; ii<numparams; ii++)
		step[ii] = 0.05;

	lastchi = feval(numparams, params, list, target, allkers);
	chi = lastchi;
	printf("starting value = %e\n", chi);

	for (ik=0; ik<30; ik++)
	{
		printf("iter %d\n", ik);
	for (ii=0; ii<numparams; ii++)
	{
		step[ii] *= 10.0;
		iter = 0;
		converged = 0;
		while (iter < 500 && !converged)
		{
			// take a step
			params[ii] += step[ii];
			// check chi
			chi = feval(numparams, params, list, target, allkers);
			if (chi > lastchi) 
			{
				step[ii] = -0.3*step[ii];
//				printf("changing step to %e: %e  %e\n", step, params[ii], chi);
			}
			if (fabs(lastchi - chi)/lastchi < tol)
			{
//				printf("reached tol after %d steps\n", iter);
				converged = 1;
			}
			iter++;
			lastchi = chi;
		}
		if (!converged)
		{
//			printf("did not converge. step=%e, iter=%d\n", step, iter);
		}
	}
	}
	printf("ending value = %e\n", chi);
}
*/
