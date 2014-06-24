#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include "/shared/cfitsio-3.24/include/fitsio.h"
#include <mpi.h>

/* Types of tags that can be parsed */
#define TAGTYPE_BOOL 0
#define TAGTYPE_INT 1
#define TAGTYPE_FLOAT 2
#define TAGTYPE_STRING 3

/* Return values for error */
#define PARSE_ERROR -1
#define TAGERR_NAME -2
#define TAGERR_TYPE -3
#define TAGERR_DATA -4

#define TASK_DONE 2
#define TASK_REQUEST 3
#define TASK_TARGET 4
#define TASK_UPDATE 5
#define TASK_ANSWER 6
#define TASK_UPDATE_DATA 7

/* used for reading command-line args */
typedef struct
{
	char *name; /* text to look for */
	int type; /* data type, see above #defines */
	void *data; /* the actual data */
} cl_tag;

/* a single kernel with included identifier */
typedef struct {
	float ***sens; /* 3d array of sensitivity for a kernel */
	int k, n; /* integer identifiers */
} kernel;

typedef struct {
	int locx, locy; // identifiers
	float flocx, flocy, size;
	int nx, ny; // size of tile on grid
	int ***ix, ***iy; // [x][y][0-3], gives grid index in x and y
	float ***cx, ***cy; // [x][y][0-15], gives coefs for each of the prev points
} interp_mask;

/* a single velocity measurement with included identifier */
typedef struct {
	float vy, ey; /* velocity and error */
	float vx, ex;
	int k, n, locx, locy; /* identifiers */
	float flocx, flocy; /* floating point location */
	float size; /* kernel size in degrees */
	kernel *ker; /* associated kernel */
	int haskernel;
	interp_mask *mask;
	int group; // 0=constant, 1=iterate
	int myindex; // in main array
} measurement;

/* kernel description, applies to all kernels */
typedef struct {
	int nx, ny, nz; // pixels in each direction
	float *z; // physical depth of each pixel in z
	float *dz;
	int maxk, maxn;
	int maxx, maxy;
	float **overlap;
	float dxdy, dx, dy;
	int numoverlap; // number of indices in measurement position
					//  that you need to move over to no longer overlap
} kernel_desc;

/* a solution point */
typedef struct {
	float vx, evx;
	float vy, evy;
	float zloc, xloc, yloc;
	float zwidth, lambda;
	float xwidth, ywidth;
	int locx, locy;
} solnpt;

typedef struct
{
	int nx, ny, nz; // number of pixels in each dimension
	float *x, *y, *z; // positions in each dimension
	float dx, dy, dxdy;
} grid;

typedef struct
{
	float x, y, z;
	float wx, wy, wz, lambda;
	float ***sens;
} target;

/* these are starting to get out of control */
typedef struct {
	float ***sens;
	float com, peak, fwhm; // some basic stats
	// more advanced stats?
	float lr; // lobe ratio. ratio of sensitivity in primary peak to non primary peaks
	float cwr; // core to wing ratio
} avgker;

/* function prototypes */

// main.c
void incremental_write(FILE *fp, solnpt *sol);
float timespec_to_sec (struct timespec *t);
void usage (char *name);

// kernels.c
void create_avgker (avgker *avg, grid *g);
void clear_avgker (avgker *avg, grid *g);
void free_avgker (avgker *avg, grid *g);
void compute_overlap(kernel_desc *allkers, int numkers, kernel **kers, int myid, int nproc);
int load_kernel_set (char *fname, char *zfname, kernel_desc *allkers, kernel*** list, int *numkers);
void free_kernel(kernel *ker);

// measurements.c
void create_grid (measurement *list, int nummeas, int numlocx, int numlocy, kernel_desc *allkers, grid *g, float resolution);
void free_grid (grid *g);
double covar (measurement *m1, measurement *m2, int group1, int group2);
void compute_answer (int numparams, float *params, measurement *mlist, measurement *****mgrid, float *vx, float *evx, float *vy, float *evy, float *norm, kernel_desc *allkers);
int link_kernels(int nummeas, int numkers, measurement *meas, kernel **kers);
int load_measurements(char *fname, measurement** list, int *nummeas, int *numlocx, int *numlocy, float **flocx, float **flocy, interp_mask **mask, int myid);
void ins_sort(int num, measurement **list, int *numunique, float **sorted, uint8_t flag);

// function.c
float feval(int numparams, float *params, measurement *list, float *target, kernel_desc *allkers);
void velocity(int numparams, float *params, measurement *list, float *vx, float *ex, float *vy, float *ey);

// optimization
void compute_cost (int numparams, float *params, measurement *mlist, target *t, avgker *avg, grid *g, kernel_desc *allkers, float *diff, float *all, float *norm);
void block_opt (int numparams, int *ind, float *params, measurement *mlist, target *t, avgker *avg, grid *g, kernel_desc *allkers, int numparams_all, measurement *****mgrid);
void inv_opt(int numparams, float *params, measurement *list, float *target, kernel_desc *allkers, float lambda);
void line_opt(int numparams, float *params, measurement *list, float *target, kernel_desc *allkers);

// interp.c
float interp_kernel (float x0);
void interp (float x, float y, float ***data, int dim, int nz, float *res);
int imin (int a, int b);
int imax (int a, int b);

// manager.c
void manager (int nproc, int numlocx, int numlocy, char *outfname, float lambda, float *flocx, float *flocy);
