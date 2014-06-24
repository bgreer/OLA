#include "header.h"

void clear_groups (measurement *mlist, int num)
{
	int ii;

	for (ii=0; ii<num; ii++)
	{
		mlist[ii].group = 0;
	}
}

double covar (measurement *m1, measurement *m2, int group1, int group2)
{
	double ret, dist, x, overlap, mode;

	overlap = 0.0;
	mode = 1.0;

	if (m1 == NULL) return 0.0;
	if (m2 == NULL) return 0.0;

	// if group is negative, ignore it
	if (group1 >= 0)
		if (m1->group != group1) return 0.0;
	if (group2 >= 0)
		if (m2->group != group2) return 0.0;

	// if different k,n
	if (m1->k != m2->k || m1->n != m2->n) return 0.0;
	if (m1->k != m2->k) return 0.0;

//	if (m1->n == m2->n) mode = 1.0;
//	else mode = (0.5 + m1->n + m2->n)/10.0;
	
	// if same position
	if (m1->locx == m2->locx && m1->locy == m2->locy) return m1->ex * m2->ex;
//	return 0.0;

	// find distance between measurement centers
	dist = sqrt(pow(m1->flocx-m2->flocx, 2.0) + pow(m1->flocy-m2->flocy, 2.0));
	if (dist > 0.5*m1->size + 0.5*m2->size) return 0.0;

	x = dist/m1->size;
//	overlap = 1 - x;
	overlap = 2.0 * (acos(x) - x*sqrt(1.0 - x*x)) / 3.14195;
	ret = m1->ex * m2->ex * overlap * mode;

	return ret;
}

// given a set of parameters, return the velocity and associated error
void compute_answer (int numparams, float *params, measurement *mlist, measurement *****mgrid, float *vx, float *evx, float *vy, float *evy, float *norm, kernel_desc *allkers)
{
	int ii, ij, ix, iy, in, ik;
	int x0, y0;
	float a, b, part;

	*vx = 0.0;
	*evx = 0.0;
	*vy = 0.0;
	*evy = 0.0;
	*norm = 0.0;
	
	for (ii=0; ii<numparams; ii++)
		*norm += params[ii];

	for (ii=0; ii<numparams; ii++)
		params[ii] /= *norm;

	for (ii=0; ii<numparams; ii++)
	{
		if (params[ii] != 0.0)
		{
			a = params[ii];
			*vx += params[ii] * mlist[ii].vx;
			*vy += params[ii] * mlist[ii].vy;
			// error w/ covariance
			for (ix=-allkers->numoverlap; ix<=allkers->numoverlap; ix++)
			{
				for (iy=-allkers->numoverlap; iy<=allkers->numoverlap; iy++)
				{
					x0 = mlist[ii].locx + ix;
					y0 = mlist[ii].locy + iy;
					if (x0 >= 0 && y0 >= 0 && x0 <= allkers->maxx && y0 <= allkers->maxy)
						for (in=0; in<=allkers->maxn; in++)
							for (ik=0; ik<=allkers->maxk; ik++)
								if (mgrid[x0][y0][in][ik] != NULL)
								{
									b = params[mgrid[x0][y0][in][ik]->myindex];
									part = a * b * covar(mlist+ii, mgrid[x0][y0][in][ik], -1, -1);
//									if (part != 0.0) printf("1 %d %f\n", ii, a * b);
									*evx += part;
								}
				}
			}
		}
	}
	*evx = sqrt(*evx);
	*evy = *evx;
}

// march through each measurement and find the appropriate kernel
// set the pointer inside the measurement to point to it
int link_kernels(int nummeas, int numkers, measurement *meas, kernel **kers)
{
	int ii, ij, counter;
	int thisn, thisk, match;

	counter = 0;
	for (ii=0; ii<nummeas; ii++)
	{
		thisn = meas[ii].n;
		thisk = meas[ii].k;

		ij = 0;
		match = -1;
		while (ij<numkers)
		{
			if (thisk == kers[ij]->k && thisn == kers[ij]->n)
			{
				// found a match
				match = ij;
				ij = numkers;
			}
			ij++;
		}
		if (match >= 0)
		{
			meas[ii].ker = kers[match];
			meas[ii].haskernel = 1;
			counter++;
		}
	}
	return counter;
}

void create_grid (measurement *list, int nummeas, int numlocx, int numlocy, kernel_desc *allkers, grid *g, float resolution)
{
	int ii, ij;
	float x0, y0, x1, y1, lts, dx, x2, x3;
	float dxmeas;

	lts = 0.0;

	// set vertical grid
	g->nz = allkers->nz;
	g->z = (float*) malloc(g->nz * sizeof(float));
	for (ii=0; ii<g->nz; ii++)
		g->z[ii] = allkers->z[ii];


	// find starting corner and largest tile size and ending corner
	for (ii=0; ii<nummeas; ii++)
	{
		if (list[ii].locx == 0) x0 = list[ii].flocx;
		if (list[ii].locx == 1) x3 = list[ii].flocx;
		if (list[ii].locy == 0) y0 = list[ii].flocy;
		if (list[ii].locx == numlocx-1) x1 = list[ii].flocx;
		if (list[ii].locy == numlocy-1) y1 = list[ii].flocy;
		if (list[ii].size > lts) lts = list[ii].size;
	}
	// shift grid start down to corner of tile
	x2 = x0;
	x0 -= lts*0.5;
	y0 -= lts*0.5;
	x1 += lts*0.5;
	y1 += lts*0.5;

	dx = allkers->dx / resolution; // in degrees
	g->dx = dx;

	// set dims
	g->nx = (x1-x0)/dx + 1;
	g->ny = (y1-y0)/dx + 1;

	// create position arrays
	g->x = (float*) malloc(g->nx * sizeof(float));
	g->y = (float*) malloc(g->ny * sizeof(float));
	for (ii=0; ii<g->nx; ii++)
		g->x[ii] = x0 + ii*dx;
	for (ii=0; ii<g->ny; ii++)
		g->y[ii] = y0 + ii*dx;

	g->dx = dx;
	g->dy = dx;
	g->dxdy = g->dx * g->dy;

	dxmeas = x3-x2;
	allkers->numoverlap = (int)(lts/dxmeas) + 1;
}

int load_measurements (char *fname, measurement** list, int *nummeas, int *numlocx, int *numlocy, float **flocx, float **flocy, interp_mask **mask)
{
	FILE *fptr;
	char instr[300];
	char *sptr;
	int ii, ij, counter, ret, numuniquex, numuniquey;
	float fread[6], *sortedx, *sortedy;
	int iread[2];
	// count number of measurements, also keep track of min/max locations
	
	fptr = fopen(fname, "r");
	// error checking
	if (fptr == NULL)
	{
		printf("Error opening measurement file: %s\n", fname);
		return -1;
	}

	// count the lines
	counter = 0;
	while(fgets(instr, 300, fptr) != NULL)
		counter++;

	// error checking?
	if (counter==0)
	{
		printf("Error: not enough measurements in %s\n", fname);
		return -2;
	}

	// make sure to return the number of measurements
	*nummeas = counter;

	// allocate space for them all
	*list = (measurement*) malloc(*nummeas * sizeof(measurement));

	// go back through and load them
	rewind(fptr);
	for (ii=0; ii<counter; ii++)
	{
		ret = fscanf(fptr,"%e\t%e\t%d\t%d\t%e\t%e\t%e\t%e\t\n", 
				fread, fread+1, iread, iread+1, fread+2, fread+3, fread+4, fread+5);
		// error checking
		if (ret != 8)
		{
			printf("Error: invalid data at line %d (%d found, 8 needed)\n", ii+1,ret);
			return -3;
		}
		// load data into list
		(*list)[ii].k = iread[0];
		(*list)[ii].n = iread[1];
		(*list)[ii].locx = -1; // not set yet
		(*list)[ii].locy = -1;
		(*list)[ii].vx = fread[2]; // zonal velocity
		(*list)[ii].ex = fread[3];
		(*list)[ii].vy = fread[4]; // meridional velocity
		(*list)[ii].ey = fread[5]; // meridional velocity error
		(*list)[ii].flocx = fread[0]; // longitude
		(*list)[ii].flocy = fread[1]; // latitude
		(*list)[ii].ker = NULL; // not set yet
		(*list)[ii].haskernel = 0; // flag for no kernel matched yet
		(*list)[ii].size = 16.0; // tile size in degreees, used to find neighbors
	}
	fclose(fptr);
	// all measurements have been loaded into memory
	
	// set integer 'loc' for each measurement by sorting floc values
	// and tagging each measurement with the sorted index
	sortedx = (float*) malloc(counter * sizeof(float));
	memset(sortedx, 0x00, counter*sizeof(float));
	ins_sort(counter, list, &numuniquex, &sortedx, 0);
	sortedy = (float*) malloc(counter * sizeof(float));
	memset(sortedy, 0x00, counter*sizeof(float));
	ins_sort(counter, list, &numuniquey, &sortedy, 1);

	printf("Found %d unique positions in x, %d in y\n", numuniquex, numuniquey);
//	for (ii=0; ii<numuniquey; ii++)
//		printf("%e\n", sortedy[ii]);
	*numlocx = numuniquex;
	*numlocy = numuniquey;

	// make space for interpolation masks
	(*mask) = (interp_mask*) malloc(numuniquex*numuniquey * sizeof(interp_mask));
	// assign identifiers
	for (ii=0; ii<numuniquex; ii++)
	{
		for (ij=0; ij<numuniquey; ij++)
		{
			(*mask)[ii*numuniquey+ij].locx = ii;
			(*mask)[ii*numuniquey+ij].locy = ij;
			(*mask)[ii*numuniquey+ij].size = 16.0;

		}
	}

	// there are now two dimensions to consider
	// we want a separate sorted list in each direction 
	// for keeping track of position

	// go through each measurement and match the floc to the sorted loc
	// assign loc to the index
	for (ii=0; ii<counter; ii++)
	{
		// search for appropriate x value
		ij = 0;
		while ((*list)[ii].flocx != sortedx[ij])
			ij++;
		(*list)[ii].locx = ij;
		// search for appropriate y value
		ij = 0;
		while ((*list)[ii].flocy != sortedy[ij])
			ij++;
		(*list)[ii].locy = ij;

		// link to an interp_mask
		ij = (*list)[ii].locx * numuniquey + (*list)[ii].locy;
		(*list)[ii].mask = &((*mask)[ij]);
	}

	// finally, every measurement has been loaded, and the 'loc' tags have been set
	
	// copy unique locations into flocx, flocy
	(*flocx) = (float*) malloc(numuniquex * sizeof(float));
	(*flocy) = (float*) malloc(numuniquey * sizeof(float));
	memcpy((*flocx), sortedx, numuniquex * sizeof(float));
	memcpy((*flocy), sortedy, numuniquey * sizeof(float));

	free(sortedx);
	free(sortedy);
	return 0;
}

// sort from lowest to highest
// insertion sort written poorly as quickly as possible from memory
// this function has been tested exactly twice.
void ins_sort(int num, measurement **list, int *numunique, float **sorted, uint8_t flag)
{
	int ii, ij, ik;
	float thisval, swap, unique;

	unique = 0;

	// flag is 0 or 1, for using flocx or flocy, respectively

	for (ij=0; ij<num; ij++)
	{
		thisval = 0.0;
		if (flag) thisval = (*list)[ij].flocy;
		else thisval = (*list)[ij].flocx;
		ii = 0;
		// find where it goes
		if (unique > 0)
			while (thisval > (*sorted)[ii] && ii < unique)
				ii++;
		// if we have hit the end, just add it
		if (ii == unique)
		{
//			printf("adding location %e at the end\n", thisval);
			(*sorted)[ii] = thisval;
			unique += 1;
		} else if (thisval != (*sorted)[ii]) {
			// if we are not at the end, and this value is unique
//			printf("inserting %e at %d\n", thisval, ii);
			// shift everything down one
			for (ik=(*numunique); ik>ii; ik--)
				(*sorted)[ik] = (*sorted)[ik-1];
			// insert
			(*sorted)[ii] = thisval;
		} else {
//			printf("non-unique value %e thrown away\n", thisval);
		}
	}
	*numunique = unique;
}

void free_grid (grid *g)
{
	free(g->x);
	free(g->y);
	free(g->z);
}
