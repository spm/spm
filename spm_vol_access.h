/* spm_vol_access.h */

/* matlab independent image access routines */

#include <sys/types.h>
#ifdef SPM_WIN32
#include "win32mmap.h"
#else
#include <sys/mman.h>
#endif

/* map mark 1 */

#define MAGIC 110494

typedef struct map
{
	double dim[3];
	double pixdim[3];
	double scale;
	int dtype;
	int off;
	caddr_t map;
	unsigned char *data;
	size_t len;
	int magic;
	int prot;
	int flags;
	pid_t pid;
}	MAP;

/* map mark 2 */

typedef struct maptype
{
	int dim[3];		/* Dimensions of the volume */
	double *scale, *offset;	/* Scalefactor and offset, such that true_intensity = vox*scale+offset */
	int dtype;		/* Data-type of volume */
	void **data;	/* Pointer to data */
	double mat[16];

	caddr_t addr;
	size_t len;
}	MAPTYPE;

int get_datasize(int type);

int resample(int m,MAPTYPE *vol,double *out,double *x, 
	double *y,double *z,int hold, double background);

int resample_d(int m,MAPTYPE *vol,double *out,
	double *gradx, double *grady, double *gradz,
	double *x, double *y, double *z,
	int hold, double background);

int slice(double *mat, double *image, int xdim1,int ydim1, 
	MAPTYPE *vol, int hold,double background);


