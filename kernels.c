#include "header.h"
/*
void characterize_kernel (kerstat *ker, kernel_desc allkers)
{
	int ii, l, r, loc;
	float ret, total, total2;

	// find center of mass
	ret = 0.0;
	total = 0.0;
	for (ii=0; ii<allkers.nz; ii++)
	{
		total += fabs(ker->avgker[ii]);
		ret += fabs(ker->avgker[ii]) * allkers.z[ii];
	}
	ker->com = ret /= total;

	// find peak
	ret = 0.0;
	total = 0.0;
	loc = -1;
	for (ii=0; ii<allkers.nz; ii++)
	{
		if (ker->avgker[ii] > total)
		{
			total = ker->avgker[ii];
			ret = allkers.z[ii];
			loc = ii;
		}
	}
	ker->peak = ret;

	// find fwhm
	l = loc+1;
	while (l > 0 && ker->avgker[l] > 0.5*total)
		l--;
	r = loc-1;
	while (r < allkers.nz-1 && ker->avgker[r] > 0.5*total)
		r++;

	ker->fwhm = allkers.z[r]-allkers.z[l];

	// find lobe ratio
	// (total sensitivity inside main lobe) / (total sensitivity outside)
	ii = loc;
	total = 0.0;
	while (ii > 0 && ker->avgker[ii] > 0.0)
	{
		total += ker->avgker[ii];
		ii--;
	}
	while (ii > 0)
	{
		total2 += fabs(ker->avgker[ii]);
		ii--;
	}
	ii = loc+1;
	while (ii < allkers.nz-1 && ker->avgker[ii] > 0.0)
	{
		total += ker->avgker[ii];
		ii++;
	}
	while (ii < allkers.nz-1)
	{
		total2 += fabs(ker->avgker[ii]);
		ii++;
	}
	ker->lr = total2/(total+total2);

	// find core to sidelobe ratio
	total = 0.0;
	for (ii=l; ii<=r; ii++)
		total += fabs(ker->avgker[ii]);
	total2 = 0.0;
	for (ii=0; ii<l; ii++)
		total2 += fabs(ker->avgker[ii]);
	for (ii=r+1; ii<allkers.nz; ii++)
		total2 += fabs(ker->avgker[ii]);
	ker->cwr = total/(total2+total);
	
}
*/

void update_avgker (avgker *avg, grid *g, int numparams, int *ind, float *params, measurement *mlist, kernel_desc *allkers)
{
	int ii, ix, iy, iz, x0, x1, y0, y1;
	float dist, par, val, x, y;
	float *res;

	res = (float*) malloc(allkers->nz * sizeof(float));

	// construct avgker
	for (ii=0; ii<numparams; ii++)
	{
		par = params[ind[ii]];
		if (par != 0.0)
		{
			// mlist[ii].ker->sens[][][];
			// coef = params[ii]
			// find bounds in x,y
			x0 = (mlist[ind[ii]].flocx - 0.5*mlist[ind[ii]].size - g->x[0])/g->dx;
			y0 = (mlist[ind[ii]].flocy - 0.5*mlist[ind[ii]].size - g->y[0])/g->dx;
			x1 = x0 + mlist[ind[ii]].size / g->dx;
			y1 = y0 + mlist[ind[ii]].size / g->dx;

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
					interp(x,y,mlist[ind[ii]].ker->sens, allkers->nx, allkers->nz, res);
					for (iz=0; iz<g->nz; iz++)
					{
						avg->sens[ix][iy][iz] += res[iz] * par;
					}
				}
			}
		}
	}
	free(res);

}

void create_avgker (avgker *avg, grid *g)
{
	int ix, iy;

	avg->sens = (float***) malloc(g->nx * sizeof(float**));
	for (ix=0; ix<g->nx; ix++)
	{
		avg->sens[ix] = (float**) malloc(g->ny * sizeof(float*));
		for (iy=0; iy<g->ny; iy++)
		{
			avg->sens[ix][iy] = (float*) malloc(g->nz * sizeof(float));
			memset(&(avg->sens[ix][iy][0]), 0x00, g->nz * sizeof(float));
		}
	}
}

void clear_avgker (avgker *avg, grid *g)
{
	int ix, iy;

	for (ix=0; ix<g->nx; ix++)
		for (iy=0; iy<g->ny; iy++)
			memset(&(avg->sens[ix][iy][0]), 0x00, g->nz * sizeof(float));

}

void free_avgker (avgker *avg, grid *g)
{
	int ix, iy;

	for (ix=0; ix<g->nx; ix++)
	{
		for (iy=0; iy<g->ny; iy++)
			free(avg->sens[ix][iy]);
		free(avg->sens[ix]);
	}
	free(avg->sens);
}

// assume the position and widths have been set
void create_target (target *t, grid *g, kernel_desc *allkers)
{
	int ix, iy, iz;
	float sum, argx, argy, argz, arg;

	sum = 0.0;

	t->sens = (float***) malloc(g->nx * sizeof(float**));
	for (ix=0; ix<g->nx; ix++)
	{
		argx = (g->x[ix] - t->x) / (t->wx*0.6);
		t->sens[ix] = (float**) malloc(g->ny * sizeof(float*));
		for (iy=0; iy<g->ny; iy++)
		{
			argy = (g->y[iy] - t->y) / (t->wy*0.6);
			t->sens[ix][iy] = (float*) malloc(g->nz * sizeof(float));
			for (iz=0; iz<g->nz; iz++)
			{
				argz = (g->z[iz] - t->z) / (t->wz*0.6);
				// gaussian
				arg = exp(-argx*argx - argy*argy - argz*argz);
				// ring
				//arg = sqrt(argx*argx + argy*argy) - 10.0;
				//arg = 1./(arg*arg+0.5+argz*argz);

				t->sens[ix][iy][iz] = arg;
				sum += arg * allkers->dz[iz] * g->dxdy;
			}
		}
	}
	// normalize
	for (ix=0; ix<g->nx; ix++)
		for (iy=0; iy<g->ny; iy++)
			for (iz=0; iz<g->nz; iz++)
				t->sens[ix][iy][iz] /= sum;

}

void free_target (target *t, grid *g)
{
	int ix, iy;

	for (ix=0; ix<g->nx; ix++)
	{
		for (iy=0; iy<g->ny; iy++)
			free(t->sens[ix][iy]);
		free(t->sens[ix]);
	}
	free(t->sens);
}

void output_target (target *t, grid *g)
{
	FILE *fp;
	int ii, ix, iy, iz;

	// output avgker
	fp = fopen("debug_t", "w");
	for (ix=0; ix<g->nx; ix++)
		for (iy=0; iy<g->ny; iy++)
			fwrite(t->sens[ix][iy], sizeof(float), g->nz, fp);
	fclose(fp);
	
}
void output_avgker (avgker *avg, grid *g)
{
	FILE *fp;
	int ii, ix, iy, iz, x0, x1, y0, y1;
	float dist;

	// output avgker
	fp = fopen("debug", "w");
	for (ix=0; ix<g->nx; ix++)
		for (iy=0; iy<g->ny; iy++)
			fwrite(avg->sens[ix][iy], sizeof(float), g->nz, fp);
	fclose(fp);
}

// this assumes the standard 3d kernel set used in ARRDI
int load_kernel_set (char *fname, char*zfname, kernel_desc *allkers, kernel*** list, int *numkers)
{
	fitsfile *fptr; // file pointer
	FILE *zfptr;
	int status = 0;
	int nx, ny, nz, nker_tot;
	long coords[4];
	float ***readin;
	float *buff, sum;
	int ii, ij, ik, counter;
	int ix, iy, iz;
	char instring [100];

	fits_open_file(&fptr, fname, READONLY, &status);
	
	// error checking
	if (status)
	{
		printf("Error opening FITS file: %s\n", fname);
		fits_report_error(stdout, status);
		return -1;
	}

	// get dimensions of kernels
	fits_read_key(fptr, TINT, "NAXIS1", &nx, NULL, &status);
	fits_read_key(fptr, TINT, "NAXIS2", &ny, NULL, &status);
	fits_read_key(fptr, TINT, "NAXIS3", &nz, NULL, &status);
	fits_read_key(fptr, TINT, "NAXIS4", &nker_tot, NULL, &status);

	// error checking
	if (status)
	{
		printf("Error reading FITS keys: %d\n", status);
		return -3;
	}

	// error checking
	if (nx < 3 || ny < 3 || nz < 3 || nker_tot < 3)
	{
		printf("Invalid kernel dimensions.\n");
		return -2;
	}

	// allocate space for the kernel list, resize later if needed
	*list = (kernel**) malloc(nker_tot * sizeof(kernel*));
	// each element is a pointer to a kernel that is floating in memory

	// copy some of this data to the kernel description
	allkers->nx = nx;
	allkers->ny = ny;
	allkers->nz = nz;
	allkers->dx = 0.9375; // dx*dy on kernel drid
	allkers->dy = allkers->dx;
	allkers->dxdy = allkers->dx * allkers->dy;
	allkers->maxk = 0;
	allkers->maxn = 0;
	allkers->maxx = 0;
	allkers->maxy = 0;
	// allocate space for the z axis
	allkers->z = (float*) malloc(nz * sizeof(float));
	allkers->dz = (float*) malloc(nz * sizeof(float));

	// temporary 3d array for reading in kernels
	readin = (float***) malloc(nx * sizeof(float**));
	for (ii=0; ii<nx; ii++)
	{
		readin[ii] = (float**) malloc(ny * sizeof(float*));
		for (ij=0; ij<ny; ij++)
			readin[ii][ij] = (float*) malloc(nz * sizeof(float));
	}
	// and a data buffer for good measure
	buff = (float*) malloc(nx * sizeof(float));

	// a note on COORDS
	// if in IDL, the array is fltarr(d1, d2, d3, d4)
	// then d1 is fastest dimension
	// a fits_read_pix will read a line [*,i,j,k]
	// with coords = {1, i, j, k}; [ uses fortran 1-indexing :( ]
	// readin is [nx][ny][nz]
	coords[0] = 1L;

	// loop through each possible kernel, check if there is anything there
	// if there is, allocate space for it, load the 2d kernel
	counter = 0;
	for (ii=0; ii<nker_tot; ii++)
	{
		// copy data from fits file to readin kernel
		coords[3] = ii+1;
		sum = 0.0;
		for (coords[2]=1L; coords[2]<=nz; coords[2]++) // loop over z
		{
			for (coords[1]=1L; coords[1]<=ny; coords[1]++) // loop over ny
			{
				fits_read_pix(fptr, TFLOAT, coords, nx, NULL, buff, NULL, &status);
				for (ij=0; ij<nx; ij++)
				{
					sum += fabs(buff[ij]); // for determining if there is a kernel here
					readin[ij][coords[1]-1][coords[2]-1] = buff[ij];
				}
			}
		}

		// if there is data in this kernel, do something with it.
		if (sum > 0.0)
		{
			(*list)[counter] = (kernel*) malloc(sizeof(kernel));
			// TODO: verify this
			(*list)[counter]->k = ii % 60;
			(*list)[counter]->n = ii/60;
			if ((*list)[counter]->n > allkers->maxn) allkers->maxn = (*list)[counter]->n;
			if ((*list)[counter]->k > allkers->maxk) allkers->maxk = (*list)[counter]->k;
			//printf("Found a kernel at %d, %d %d\n", ii, (*list)[counter]->k, (*list)[counter]->n);
			// allocate space in the kernel
			(*list)[counter]->sens = (float***) malloc(nx * sizeof(float**));
			for (ij=0; ij<nx; ij++)
			{
				(*list)[counter]->sens[ij] = (float**) malloc(ny * sizeof(float*));
				for (ik=0; ik<ny; ik++)
					(*list)[counter]->sens[ij][ik] = (float*) malloc(nz * sizeof(float));
			}
			// copy readin data into kernel sensitivity
			for (ix=0; ix<nx; ix++)
				for (iy=0; iy<ny; iy++)
					for (iz=0; iz<nz; iz++)
						(*list)[counter]->sens[ix][iy][iz] = readin[ix][iy][iz];
			counter++;
		}
	}
	// close FITS file
	fits_close_file(fptr, &status);
	*numkers = counter;

	// go ahead and read in the physical z levels
	zfptr = fopen(zfname, "r");
	if (zfptr == NULL)
	{
		printf("Error opening file %s\n", zfname);
		return -4;
	}
	for (ii=0; ii<nz; ii++)
	{
		fgets(instring, 100, zfptr);
		// split string
		ij = 0;
		while (instring[ij] == ' ') ij++;
		while (instring[ij] != ' ') ij++;
		instring[ij] = '\0'; // atof uses null char as stopping point
		allkers->z[ii] = atof(instring);
		allkers->dz[ii] = atof(instring+ij+1);
	}

	fclose(zfptr);


	// free temporary memory
	free(buff);
	for (ii=0; ii<nx; ii++)
	{
		for (ij=0; ij<ny; ij++)
			free(readin[ii][ij]);
		free(readin[ii]);
	}
	free(readin);


	return 0;
}

// precompute the overlap between kernels located at the same position
// want to address an array overlap[ind1][ind2]
// where ind1 = n1*maxk + k1
// since this is pre-collapsed, compute on kernel grid, not solution grid
void compute_overlap(kernel_desc *allkers, int numkers, kernel **kers, int myid, int nproc)
{
	int nummodes;
	int ii, ij, ind1, ind2;
	int ix, iy, iz;
	float dist;
	// mpi:
	int start, end;

	nummodes = (allkers->maxn+1) * (allkers->maxk+1);

	/* example:
	 * k1 = 25
	 * n1 = 0
	 * k2 = 10
	 * n2 = 0
	 * result = 0.003816
	 * current = 0.000890
	 */ 

	// allocate space for overlap matrix
	allkers->overlap = (float**) malloc(nummodes * sizeof(float*));
	for (ii=0; ii<nummodes; ii++)
		allkers->overlap[ii] = (float*) malloc(nummodes * sizeof(float));

	// divide task
	start = myid*numkers/nproc;
	end = (myid+1)*numkers/nproc - 1;

	for (ii=start; ii<=end; ii++)
	{
		ind1 = kers[ii]->n*allkers->maxk + kers[ii]->k;
		for (ij=0; ij<numkers; ij++)
		{
			ind2 = kers[ij]->n*allkers->maxk + kers[ij]->k;
			allkers->overlap[ind1][ind2] = 0.0;
			// sum kers[ii] * kers[ij] * dz
			for (ix=0; ix<allkers->nx; ix++)
			{
				for (iy=0; iy<allkers->ny; iy++)
				{
					for (iz=0; iz<allkers->nz; iz++)
						allkers->overlap[ind1][ind2] += 
							kers[ii]->sens[ix][iy][iz] * kers[ij]->sens[ix][iy][iz] * 
							allkers->dz[iz] * allkers->dxdy;
				}
			}
			/*
			if (kers[ii]->k == 25 && kers[ii]->n == 0)
			{
				printf("%d %d %f\n", kers[ij]->k, kers[ij]->n, );
				exit(0);
			}*/
		}
	}

	// share answers
	for (ii=0; ii<nproc; ii++)
	{
		start = ii*numkers/nproc;
		end = (ii+1)*numkers/nproc - 1;
		for (ij=start; ij<=end; ij++)
		{
			ind1 = kers[ij]->n*allkers->maxk + kers[ij]->k;
			MPI_Bcast(allkers->overlap[ind1], nummodes, MPI_FLOAT, ii, MPI_COMM_WORLD);
		}
	}

}

// TODO: this, later
void free_kernel(kernel *ker)
{
	
}
