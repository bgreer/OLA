#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "header.h"

/* General flow:
 * 1 - Load measurements
 * 2 - Load associated kernels
 * 3 - Link measurements to kernels
 * 4 - Create solution grid, each point solved independently
 * For each grid point,
 * 	 5 - Determine what measurements to include
 * 	 6 - Create inversion grid (spans included neighbors)
 * 	 7 - Create KTK matrix by placing measurement kernels on inversion grid
 * 	 8 - Invert matrix, produce coefficients
 * 	 9 - Apply 
 *
 * trying to minimize the quantity:
 * 		sum_xyz [ T(r) - sum_i a_i K_i ]^2 + lambda sum_i [a_i sigma_i]^2
 * do this by taking derivatives wrt each a_j, setting all to zero, solving together
 * unfortunately, this cannot be done in one blow if there are too many a's
 * so we take the derivative wrt some subset of a_j, solve that block together
 * we are still solving as a = (K^T K + R)^-1 K^T T
 * but with a smaller matrix to invert each iteration
 * we still need to compute K^T K + R all together instead of K, then computing rest
 * the blocksize is the number of parameters to include in each iteration
 * there might be some optimal value, something like 250.
 * also, figure out best blocking algorithm
 *
 * the kernel matrix is nparams by ngridpts. need to know number of grid points
 * go with grid that kernels are computed on, usually 125 pixels in depth, 17 in x and y
 * PROBLEM: arbitrary tile overlap means crazy grid. need to interpolate kernels
 *  center of 17x17 kernel is at pixel 8,8 (0-indexed)
 * need to determine optimal grid based on tile overlap.
 * use tile centers as initial grid. probably want to at least resolve the kernel grid
 * what if tile center does not line up on adjacent kernel grid?
 * start grid with corner tile, use kernel grid from that everywhere
 * interpolate all other kernels on to this grid
*/

int main (int argc, char *argv[])
{
	measurement *meas, *mlist;
	measurement *****mgrid; // mgrid[x][y][n][k] is a pointer to a measurement
	grid g;
	kernel_desc allkers;
	kernel **kers;
	target t;
	avgker avg;
	interp_mask *mask;
	cl_tag tag;
	// non-structs:
	int nummeas, numlocx, numlocy;
	int numkers;
	float lambda, *flocx, *flocy, x, y;
	int ii, ij, ik, il, niter, iter, ix, iy, iz;
	char infname[256], outfname[256];
	struct timespec time_start, time_end;
	float dt, resolution;
	// iteration:
	int numparams_all, blocksize, maxblocksize, *ind;
	int listsize;
	int *listx, *listy;
	float *listval, **diff, sum, dist;
	float *params, *params_diff, *params_old;
	// answer output:
	FILE *fp;
	float vx, vy, evx, evy, dev, cost, norm, norm2;
	// blocking:
	int thisx, thisy, redo;
	// mpi:
	int nproc, myid, newtask;
	int mpisend;
	float mpifltbuffer[4];
	MPI_Status stat;

	// set up MPI stuff
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	printf("Processor %d of %d checking in.\n", myid, nproc);
	if (nproc == 1)
	{
		printf("ERROR: Need at least 2 MPI processes.\n");
		exit(-1);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	// use some command-line arguments
	tag.name = "-in";
	tag.type = TAGTYPE_STRING;
	tag.data = infname;
	if (!parse_params(argc, argv, 1, &tag))
		{usage(argv[0]);return -1;}
	tag.name = "-out";
	tag.type = TAGTYPE_STRING;
	tag.data = outfname;
	if (!parse_params(argc, argv, 1, &tag))
		{usage(argv[0]);return -1;}
	lambda = 6.0;
	tag.name = "-reg";
	tag.type = TAGTYPE_FLOAT;
	tag.data = &lambda;
	parse_params(argc, argv, 1, &tag);
	resolution = 2.0;
	tag.name = "-res";
	tag.type = TAGTYPE_FLOAT;
	tag.data = &resolution;
	parse_params(argc, argv, 1, &tag);

	// print some run details
	if (myid == 0)
	{
		printf("\tUsing input file %s\n", infname);
		printf("\tUsing output file %s\n", outfname);
		printf("\tUsing lambda = %f\n", lambda);
		printf("\tGrid:Kernel resolution = %f\n", resolution);
	}


	// load measurements from file (parallel)
	load_measurements(infname, &meas, &nummeas, &numlocx, &numlocy,
			&flocx, &flocy, &mask, myid);


	// load kernels
	load_kernel_set("idl/kers_128.fits", "idl/z_txt", &allkers, &kers, &numkers);
	if (myid == 0)
	{
		printf("Number of kernels loaded: %d\n", numkers);
		printf("Kernel dims: %d %d %d\n", allkers.nx, allkers.ny, allkers.nz);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	
	
	// pre-compute kernel overlap
	if (myid == 0) printf("Precomputing 1-D Overlap Matrix..\n");
	compute_overlap(&allkers, numkers, kers, myid, nproc); // parallel

	// set up grid that kernels will be placed on
	create_grid(meas, nummeas, numlocx, numlocy, &allkers, &g, resolution);

	// pre-compute interpolation
	for (ii=0; ii<numlocx; ii++)
		for (ij=0; ij<numlocy; ij++)
			precompute_interp(mask+(ii*numlocx)+ij);

	// link each measurement to a kernel
	numparams_all = link_kernels(nummeas, numkers, meas, kers);
	if (myid == 0) printf("Number of measurements linked to kernels: %d\n", numparams_all);

	// copy linked measurements into new list
	mlist = (measurement*) malloc(numparams_all * sizeof(measurement));
	ij = 0;
	for (ii=0; ii<numparams_all; ii++)
	{
		while (!meas[ij].haskernel) ij++;
		// copy meas[ij] into mlist[ii]
		memcpy(mlist+ii, meas+ij, sizeof(measurement));
		ij++;
	}
	free(meas);
	allkers.maxx = numlocx-1;
	allkers.maxy = numlocy-1;

	// create measurement grid
	mgrid = (measurement*****) malloc(numlocx * sizeof(measurement****));
	for (ix=0; ix<numlocx; ix++)
	{
		mgrid[ix] = (measurement****) malloc(numlocy * sizeof(measurement***));
		for (iy=0; iy<numlocy; iy++)
		{
			mgrid[ix][iy] = (measurement***) malloc((allkers.maxn+1) * sizeof(measurement**));
			for (ii=0; ii<=allkers.maxn; ii++)
			{
				mgrid[ix][iy][ii] = (measurement**) malloc((allkers.maxk+1) * sizeof(measurement*));
				for (ij=0; ij<=allkers.maxk; ij++)
					mgrid[ix][iy][ii][ij] = NULL;
			}
		}
	}
	// load measurements into grid
	if (myid == 0) printf("Loading Measurement Grid..\n");
	for (ii=0; ii<numparams_all; ii++)
	{
		mgrid[mlist[ii].locx][mlist[ii].locy][mlist[ii].n][mlist[ii].k] = mlist+ii;
		mlist[ii].myindex = ii;
	}


	// create averaging kernel
	if (myid == 0) printf("Initializing Averaging Kernel..\n");
	create_avgker(&avg, &g);
	
	// set up arrays for iterating inversion
	params = (float*) malloc(numparams_all * sizeof(float));
	params_diff = (float*) malloc(numparams_all * sizeof(float));
	params_old = (float*) malloc(numparams_all * sizeof(float));
	maxblocksize = 3000;
	if (myid == 0) printf("Using an max iteration blocksize of %d.\n", maxblocksize);
	ind = (int*) malloc(maxblocksize * sizeof(int));


	// initialize parameters
	listsize = 10;
	listx = (int*) malloc(listsize * sizeof(int));
	listy = (int*) malloc(listsize * sizeof(int));
	listval = (float*) malloc(listsize * sizeof(float));
	diff = (float**) malloc(g.nx * sizeof(float*));
	for (ix=0; ix<g.nx; ix++)
		diff[ix] = (float*) malloc(g.ny * sizeof(float));


	dt = 0.0;

	// split off rank 0 for managing
	if (myid != 0)
	{
		printf("P%d\tRequesting task..\n", myid);
		// request first task
		MPI_Send(&mpisend, 1, MPI_INT, 0, TASK_REQUEST, MPI_COMM_WORLD);
		printf("P%d\tWaiting for response..\n", myid);
		// get answer to request
		MPI_Recv(&mpisend, 1, MPI_INT, 0, TASK_REQUEST, MPI_COMM_WORLD, &stat);
		if (mpisend == 1)
		{
			newtask = 1;
		} else {
			newtask = 0;
		}
		
		// while there are things to do,
		while (newtask)
		{
		// get target from manager
		MPI_Recv(&t, sizeof(target), MPI_BYTE, 0, TASK_TARGET, MPI_COMM_WORLD, &stat);
		printf("P%d\tCreating target at %f %f %f\n", myid, t.x, t.y, t.z);
		create_target(&t, &g, &allkers);

	//	output_target(&t, &g);
		memset(params, 0x00, numparams_all * sizeof(float));
		memset(params_diff, 0x00, numparams_all * sizeof(float));
		memset(params_old, 0x00, numparams_all * sizeof(float));

		
		thisx = 0;
		thisy = 0;
		niter = 5;
		// START ITERATION
		for (iter=0; iter<niter; iter++)
		{
			printf("P%d\tIteration %d at (%f, %f, %f).\n", myid, iter+1, t.x, t.y, t.z);
	
			// what positions would be best to iterate?
			// loop through each position, determine how 'bad' the profile is there
			// sort positions by badness
	
			// initialize badness
			memset(listval, 0x00, listsize * sizeof(float));
			// compute difference
			for (ix=0; ix<g.nx; ix++)
			{
				for (iy=0; iy<g.ny; iy++)
				{
					diff[ix][iy] = 0.0;
					for (iz=0; iz<g.nz; iz++)
						diff[ix][iy] += t.sens[ix][iy][iz] - avg.sens[ix][iy][iz];
				}
			}
			// loop through positions
			for (ii=0; ii<numlocx; ii++)
			{
				for (ij=0; ij<numlocy; ij++)
				{
					sum = 0.0;
					for (ix=0; ix<g.nx; ix++)
					{
						for (iy=0; iy<g.ny; iy++)
						{
							dist = sqrt(pow(g.x[ix]-flocx[ii],2.0)+pow(g.y[iy]-flocy[ij],2.0));
							if (dist < 8.0)
								sum += fabs(diff[ix][iy]);
						}
					}
					// insertion sort into list
					ik = 0;
					while (ik < listsize)
					{
						if (listval[ik] >= sum) // move along
							ik++;
						else // insert
						{
							// shift everyone else down one
							for (il=listsize-1; il>ik; il--)
							{
								listval[il] = listval[il-1];
								listx[il] = listx[il-1];
								listy[il] = listy[il-1];
							}
							// insert at ik
							listval[ik] = sum;
							listx[ik] = ii;
							listy[ik] = ij;
							ik = listsize;
						}
					}
				}
			}
	
			// do most important locations first
			for (ii=0; ii<((iter*2+8<listsize)?iter*2+8:listsize); ii++)
			{
				clear_groups(mlist, numparams_all);
				blocksize = 0;
				ij = 0;
				thisx = listx[ii];// + rand()%7 - 3;
				thisy = listy[ii];// + rand()%7 - 3;
/*				if (rand()%10 == 0)
				{
					thisx += rand()%7-3;
					thisy += rand()%7-3;
				}*/
				while (blocksize < maxblocksize && ij < numparams_all)
				{
					// TODO: optimize this using mgrid[][][][]
					if (mlist[ij].locx == thisx && mlist[ij].locy == thisy)
					{
						ind[blocksize] = ij;
						mlist[ij].group = 1;
						blocksize++;
					}
					ij++;
				}
				if (blocksize == 0) continue;
	//			printf("P%d\tI%d: %d %d %f\n", myid, ii, thisx, thisy, listval[ii]);
				memcpy(params_old, params, numparams_all * sizeof(float));
				// iterate block
				clock_gettime(CLOCK_MONOTONIC, &time_start);
				block_opt(blocksize, ind, params, mlist, &t, &avg, &g, &allkers, numparams_all, mgrid);
				// update averaging kernel
//				printf("P%d\tUpdating Averaging Kernel..\n", myid);
				// computing param difference
				for (ij=0; ij<numparams_all; ij++)
					params_diff[ij] = params[ij] - params_old[ij];
				update_avgker(&avg, &g, blocksize, ind, params_diff, mlist, &allkers);
				clock_gettime(CLOCK_MONOTONIC, &time_end);
				dt += (timespec_to_sec(&time_end) - timespec_to_sec(&time_start));
		/*	
				output_avgker(&avg, &g);
				// send update to manager
				compute_answer(numparams_all, params, mlist, mgrid, &vx, &evx, &vy, &evy, &norm, &allkers);
				compute_cost(numparams_all, params, mlist, &t, &avg, &g, &allkers, &dev, &cost, &norm2);
				MPI_Send(&mpisend, 1, MPI_INT, 0, TASK_UPDATE, MPI_COMM_WORLD);
				mpifltbuffer[0] = vx;
				mpifltbuffer[1] = evx;
				mpifltbuffer[2] = dev;
				MPI_Send(mpifltbuffer, 3, MPI_FLOAT, 0, TASK_UPDATE_DATA, MPI_COMM_WORLD);
				*/
			}
		}
		// END ITERATION
		printf("P%d\tTotal Inversion Time: %8.2fs\n", myid, dt);
		printf("P%d\tComputing Answer..\n", myid);
		vx = 0.0; vy = 0.0; evx = 0.0; evy = 0.0;
		compute_answer(numparams_all, params, mlist, mgrid, &vx, &evx, &vy, &evy, &norm, &allkers);

		// free target function
		free_target(&t, &g);
		clear_avgker(&avg, &g);
		// tell manager we are done
		MPI_Send(&mpisend, 1, MPI_INT, 0, TASK_DONE, MPI_COMM_WORLD);
		// package and send answer
		mpifltbuffer[0] = vx;
		mpifltbuffer[1] = evx;
		mpifltbuffer[2] = vy;
		mpifltbuffer[3] = evy;
		MPI_Send(mpifltbuffer, 4, MPI_FLOAT, 0, TASK_ANSWER, MPI_COMM_WORLD);
		// request new task
		MPI_Send(&mpisend, 1, MPI_INT, 0, TASK_REQUEST, MPI_COMM_WORLD);
		// get answer to request
		MPI_Recv(&mpisend, 1, MPI_INT, 0, TASK_REQUEST, MPI_COMM_WORLD, &stat);
		if (mpisend == 1)
		{
			newtask = 1;
		} else {
			newtask = 0;
		}
		}


	} else {
		manager(nproc, numlocx, numlocy, outfname, lambda, flocx, flocy);
	}

	// free memory
	for (ii=0; ii<numkers; ii++)
	{
		for (ix=0; ix<allkers.nx; ix++)
		{
			for (iy=0; iy<allkers.ny; iy++)
				free(kers[ii]->sens[ix][iy]);
			free(kers[ii]->sens[ix]);
		}
		free(kers[ii]);
	}
	free(kers);
	free_avgker(&avg, &g);
	free(allkers.z);
	free(allkers.dz);
	for (ii=0; ii<(allkers.maxn+1) * (allkers.maxk+1); ii++)
		free(allkers.overlap[ii]);
	free(allkers.overlap);
	free_grid(&g);
	free(flocx);
	free(flocy);
	free(ind);
	free(params);
	free(params_diff);
	free(mlist);
	free(listx);
	free(listy);
	free(listval);
	for (ix=0; ix<g.nx; ix++)
		free(diff[ix]);
	free(diff);
	// free mgrid
	for (ix=0; ix<numlocx; ix++)
	{
		for (iy=0; iy<numlocy; iy++)
		{
			for (ii=0; ii<=allkers.maxn; ii++)
				free(mgrid[ix][iy][ii]);
			free(mgrid[ix][iy]);
		}
		free(mgrid[ix]);
	}
	free(mgrid);


	printf("P%d\t Ready to exit.\n", myid);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return EXIT_SUCCESS;
}

void incremental_write(FILE *fp, solnpt *sol)
{
	fprintf(fp, "%e\t%e\t%e\t%e\t%e\t%e\t%e\n", sol->xloc, sol->yloc, sol->zloc, sol->vx, sol->evx, sol->vy, sol->evy);
	fflush(fp);
}

float timespec_to_sec (struct timespec *t)
{
	return (float)(t->tv_sec) + 1e-9 * (uint32_t)(t->tv_nsec);
}

void usage (char *name)
{
	printf("Usage:\n\t%s\n", name);
}
