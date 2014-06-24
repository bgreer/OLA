
#include "header.h"

void manager (int nproc, int numlocx, int numlocy, char *outfname, float lambda, float *flocx, float *flocy)
{
	FILE *fp, *fp2;
	int ix, iy, iz, ii;
	int *currenttask;
	int numtasks, task, procsrunning, flag;
	solnpt *sol;
	MPI_Status stat;
	float mpifltbuffer[4];
	target t;
	int mpisend;

	printf("MANAGER START\n");
	// for keeping track of what each proc is doing
	currenttask = (int*) malloc(nproc * sizeof(int));
	memset(currenttask, 0x00, nproc * sizeof(int));

	// create a grid of solution points
	numtasks = numlocx * numlocy * 1;
	sol = (solnpt*) malloc(numtasks * sizeof(solnpt));
	task = 0;
	for (iz=0; iz<1; iz++)
	{
		for (iy=0; iy<numlocy; iy++)
		{
			for (ix=0; ix<numlocx; ix++)
			{
				sol[task].xloc = flocx[ix];
				sol[task].yloc = flocy[iy];
				sol[task].zloc = 30.0;// - iz*4.0;
				sol[task].xwidth = 0.5;
				sol[task].ywidth = sol[task].xwidth;
				sol[task].zwidth = 1.0 + sol[task].zloc*0.15;
				sol[task].lambda = lambda + 0.12*sol[task].zloc;
				sol[task].locx = ix;
				sol[task].locy = iy;
				task++;
			}
		}
	}
	task = 0;

	printf("MANAGER starting %d tasks.\n", numtasks);

	fp = fopen(outfname, "w");
	fp2 = fopen("converge", "w");

	procsrunning = nproc-1;
	// while there are solution points to do, dish them out to processors
	while(procsrunning>0)
	{
		// loop through processors, checking for messages
		for (ii=1; ii<nproc; ii++)
		{
			// check for task update
			flag = 0;
			MPI_Iprobe(ii, TASK_UPDATE, MPI_COMM_WORLD, &flag, &stat);
			if (flag)
			{
				MPI_Recv(&mpisend, 1, MPI_INT, ii, TASK_UPDATE, MPI_COMM_WORLD, &stat);
				MPI_Recv(mpifltbuffer, 3, MPI_FLOAT, ii, TASK_UPDATE_DATA, MPI_COMM_WORLD, &stat);
				printf("MANAGER received update from proc %d (task %d)\n", ii, currenttask[ii]);
				fprintf(fp2, "%f\t%f\t%f\n", mpifltbuffer[0], mpifltbuffer[1], mpifltbuffer[2]);
				fflush(fp2);
			}
			// check for task completed
			flag = 0;
			MPI_Iprobe(ii, TASK_DONE, MPI_COMM_WORLD, &flag, &stat);
			if (flag)
			{
				MPI_Recv(&mpisend, 1, MPI_INT, ii, TASK_DONE, MPI_COMM_WORLD, &stat);
				// receive some data...
				MPI_Recv(mpifltbuffer, 4, MPI_FLOAT, ii, TASK_ANSWER, MPI_COMM_WORLD, &stat);
				// unpack data
				sol[currenttask[ii]].vx = mpifltbuffer[0];
				sol[currenttask[ii]].evx = mpifltbuffer[1];
				sol[currenttask[ii]].vy = mpifltbuffer[2];
				sol[currenttask[ii]].evy = mpifltbuffer[3];
				incremental_write(fp, &(sol[currenttask[ii]]));
				currenttask[ii] = -1;
			}
			// check for task request
			flag = 0;
			MPI_Iprobe(ii, TASK_REQUEST, MPI_COMM_WORLD, &flag, &stat);
			if (flag)
			{
				// receiving message to get rid of it
				MPI_Recv(&mpisend, 1, MPI_INT, ii, TASK_REQUEST, MPI_COMM_WORLD, &stat);
				if (task<numtasks)
				{
					// tell if there is a new task to do
					mpisend = 1;
					MPI_Send(&mpisend, 1, MPI_INT, ii, TASK_REQUEST, MPI_COMM_WORLD);
					// send that task..
					printf("MANAGER received request from %d for task, sending task %d\n", ii, task);
					t.x = sol[task].xloc;
					t.y = sol[task].yloc;
					t.z = sol[task].zloc;
					t.wx = sol[task].xwidth;
					t.wy = sol[task].ywidth;
					t.wz = sol[task].zwidth;
					t.lambda = sol[task].lambda;
					MPI_Send(&t, sizeof(target), MPI_BYTE, ii, TASK_TARGET, MPI_COMM_WORLD);
					currenttask[ii] = task;
					task++;
				} else { // nothing else to do, send no task
					printf("MANAGER received request from %d for task, denying..\n", ii);
					mpisend = 0;
					MPI_Send(&mpisend, 1, MPI_INT, ii, TASK_REQUEST, MPI_COMM_WORLD);
					currenttask[ii] = -1;
					procsrunning--;
				}
			}
		}
	}
	fclose(fp2);
	fclose(fp);
	free(sol);

	printf("MANAGER DONE\n");

	free(currenttask);
}
