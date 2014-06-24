#include "/home8/begr7169/SOFTWARE/lapack-3.4.0/lapacke/include/lapacke.h"

int main (int argc, char *argv[])
{
	double *ktk, *rscale, *cscale, *work, *ans, *ans2;
	int *piv;
	double row, col, amax, temp;
	lapack_int n, lda, info, worksize;
	int ii, ij, iz, numparams;

	numparams = atoi(argv[1]);

	printf("N = %d\n", numparams);
	
	ktk = (double*) malloc(numparams * numparams * sizeof(double));
	piv = malloc(numparams*numparams*sizeof(int));
	rscale = malloc(numparams*sizeof(double));
	cscale = malloc(numparams*sizeof(double));

	// load ktk
	printf("Loading matrix..\n");
	for (ii=0; ii<numparams; ii++)
	{
		for (ij=0; ij<numparams; ij++)
		{
			ktk[ii*numparams+ij] = (rand()%1000)/1000.0;
		}
	}

	// pre-conditioning
	printf("Preconditioning..\n");
	n = numparams;
	lda = numparams;
	dgeequ_(&n, &n, ktk, &lda, rscale, cscale, &row, &col, &amax, &info);
	free(rscale);
	free(cscale);

	// lu decomposition
	printf("LU decomposition..\n");
	n = numparams;
	lda = numparams;
	dgetrf_(&n, &n, ktk, &lda, piv, &info);
	if (info) printf("\tLU decomp status: %d\n", info);

	// inversion
	printf("Inversion..\n");
	worksize = 32*n;
	work = malloc(worksize*sizeof(double));
	dgetri_(&n, ktk, &lda, piv, work, &worksize, &info);
	if (info) printf("\tInversion status: %d\n", info);
	printf("Done.\n");

	free(work);
	free(piv);
	free(ktk);
}

