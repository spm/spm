#ifndef lint
static char sccsid[]="%W% John Ashburner %E%";
#endif

/* matlab dependent high level data access and map manipulation routines */

#include <math.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include "spm_sys_deps.h"
#include "spm_mapping.h"
#include "spm_datatypes.h"

/**************************************************************************/

void free_maps(MAPTYPE *maps, int n)
{
	int j;
	for(j=0; j<n; j++)
	{
		if (maps[j].addr)
		{
			#ifdef SPM_WIN32
			(void)unmap_file(maps[j].addr);
			#else
			(void)munmap(maps[j].addr, maps[j].len);
			#endif
			maps[j].addr=0;
		}
		if (maps[j].data)
		{
			(void)mxFree((char *)maps[j].data);
			maps[j].data=0;
		}
		if (maps[j].scale)
		{
			(void)mxFree((char *)maps[j].scale);
			maps[j].scale=0;
		}
		if (maps[j].offset)
		{
			(void)mxFree((char *)maps[j].offset);
			maps[j].offset=0;
		}
	}
	(void)mxFree((char *)maps);
}

/**************************************************************************/

static void get_map_dat(int i, const mxArray *ptr, MAPTYPE *maps)
{
	mxArray *tmp;
	double *pr;
	int num_dims, j, t, dtype;
	const int *dims;
	unsigned char *dptr;

	tmp=mxGetField(ptr,i,"dat");
	if (tmp == (mxArray *)0)
	{
		free_maps(maps,i);
		mexErrMsgTxt("Cant find dat.");
	}
	if      (mxIsDouble(tmp)) dtype = SPM_DOUBLE;
	else if (mxIsSingle(tmp)) dtype = SPM_FLOAT;
	else if (mxIsInt32 (tmp)) dtype = SPM_SIGNED_INT;
	else if (mxIsUint32(tmp)) dtype = SPM_UNSIGNED_INT;
	else if (mxIsInt16 (tmp)) dtype = SPM_SIGNED_SHORT;
	else if (mxIsUint16(tmp)) dtype = SPM_UNSIGNED_SHORT;
	else if (mxIsInt8  (tmp)) dtype = SPM_SIGNED_CHAR;
	else if (mxIsUint8 (tmp)) dtype = SPM_UNSIGNED_CHAR;
	else
	{
		free_maps(maps,i);
		mexErrMsgTxt("Unknown volume datatype.");
	}

	num_dims = mxGetNumberOfDimensions(tmp);
	if (num_dims > 3)
	{
		free_maps(maps,i);
		mexErrMsgTxt("Too many dimensions.");
	}

	dims     = mxGetDimensions(tmp);
	for(j=0; j<num_dims; j++)
		maps[i].dim[j]=dims[j];
	for(j=num_dims; j<3; j++)
		maps[i].dim[j]=1;

	tmp=mxGetField(ptr,i,"dim");
	if (tmp != (mxArray *)0)
	{
		if (mxGetM(tmp)*mxGetN(tmp) != 4)
		{
			free_maps(maps,i);
			mexErrMsgTxt("Wrong sized dim.");
		}
		pr = mxGetPr(tmp);
		if (maps[i].dim[0] != (int)fabs(pr[0]) ||
		    maps[i].dim[1] != (int)fabs(pr[1]) ||
		    maps[i].dim[2] != (int)fabs(pr[2]))
		{
			free_maps(maps,i);
			mexErrMsgTxt("Incompatible volume dimensions in dim.");
		}
		if (maps[i].dtype != (int)fabs(pr[3]))
		{
			free_maps(maps,i);
			mexErrMsgTxt("Incompatible datatype in dim.");
		}
	}

	maps[i].addr      = 0;
	maps[i].len       = 0;
	maps[i].dtype  = dtype;
	maps[i].data   = (void  **)mxCalloc(maps[i].dim[2],sizeof(void *));
	maps[i].scale  = (double *)mxCalloc(maps[i].dim[2],sizeof(double));
	maps[i].offset = (double *)mxCalloc(maps[i].dim[2],sizeof(double));


	t     = maps[i].dim[0]*maps[i].dim[1]*get_datasize(maps[i].dtype)/8;
	dptr  = (unsigned char *)mxGetPr(tmp);

	tmp=mxGetField(ptr,i,"pinfo");
	if (tmp != (mxArray *)0)
	{
		if (mxGetM(tmp) != 2 || (mxGetN(tmp) != 1 && mxGetN(tmp) != maps[i].dim[2]))
		{
			free_maps(maps,i+1);
			mexErrMsgTxt("Wrong sized pinfo.");
		}
		pr = mxGetPr(tmp);
		if (mxGetN(tmp) == 1)
			for(j=0; j<maps[i].dim[2]; j++)
			{
				maps[i].scale[j]  = pr[0];
				maps[i].offset[j] = pr[1];
				maps[i].data[j]   = &(dptr[j*t]);
			}
		else
			for(j=0; j<maps[i].dim[2]; j++)
			{
				maps[i].scale[j]  = pr[0+j*2];
				maps[i].offset[j] = pr[1+j*2];
				maps[i].data[j]   = &(dptr[j*t]);
			}
	}
	else
		for(j=0; j<maps[i].dim[2]; j++)
		{
			maps[i].scale[j]  = 1.0;
			maps[i].offset[j] = 0.0;
			maps[i].data[j]   = &(dptr[j*t]);
		}


	tmp=mxGetField(ptr,i,"mat");
	if (tmp != (mxArray *)0)
	{
		if (mxGetM(tmp) != 4 || mxGetN(tmp) != 4)
		{
			free_maps(maps,i+1);
			mexErrMsgTxt("Wrong sized mat.");
		}
		pr = mxGetPr(tmp);
		for(j=0; j<16; j++)
			maps[i].mat[j] = pr[j];
	}
	else
	{
		for(j=0; j<16; j++)
			maps[i].mat[j] = 0.0;
		for(j=0; j<4; j++)
			maps[i].mat[j + j*4] = 1.0;
	}
}

/**************************************************************************/

static void get_map_file(int i, const mxArray *ptr, MAPTYPE *maps)
{
	int j;
	mxArray *tmp;
	double *pr;
	int dsize, off;

	tmp=mxGetField(ptr,i,"dim");
	if (tmp == (mxArray *)0)
	{
		free_maps(maps,i);
		mexErrMsgTxt("Cant find dim.");
	}
	if (mxGetM(tmp)*mxGetN(tmp) != 4)
	{
		free_maps(maps,i);
		mexErrMsgTxt("Wrong sized dim.");
	}
	pr = mxGetPr(tmp);
	maps[i].dim[0] = (int)fabs(pr[0]);
	maps[i].dim[1] = (int)fabs(pr[1]);
	maps[i].dim[2] = (int)fabs(pr[2]);
	maps[i].dtype  = (int)fabs(pr[3]);
	maps[i].data   = (void  **)mxCalloc(maps[i].dim[2],sizeof(void *));
	maps[i].scale  = (double *)mxCalloc(maps[i].dim[2],sizeof(double));
	maps[i].offset = (double *)mxCalloc(maps[i].dim[2],sizeof(double));

	dsize = get_datasize(maps[i].dtype);
	if (dsize==0)
	{
		free_maps(maps,i+1);
		mexErrMsgTxt("Unknown datatype.");
	}


	tmp=mxGetField(ptr,i,"fname");
	if (tmp == (mxArray *)0)
	{
		free_maps(maps,i+1);
		mexErrMsgTxt("Cant find fname.");
	}
	if (mxIsChar(tmp))
	{
		int buflen;
		char *buf;
		int fd;
		struct stat stbuf;
		buflen = mxGetN(tmp)*mxGetM(tmp)+1;
		buf = mxCalloc(buflen,sizeof(char));
		if (mxGetString(tmp,buf,buflen))
		{
			mxFree(buf);
			free_maps(maps,i+1);
			mexErrMsgTxt("Cant get filename.");
		}
		if ((fd = open(buf, O_RDONLY)) == -1)
		{
			mxFree(buf);
			free_maps(maps,i+1);
			mexErrMsgTxt("Cant open image file.");
		}
		if (fstat(fd, &stbuf) == -1)
		{
			(void)close(fd);
			mxFree(buf);
			free_maps(maps,i+1);
			mexErrMsgTxt("Cant stat image file.");
		}
		maps[i].len = stbuf.st_size;
		#ifdef SPM_WIN32
		(void)close(fd);
		maps[i].addr = map_file(buf, (caddr_t)0, maps[i].len,
			PROT_READ, MAP_SHARED, (off_t)0);
		#else
		maps[i].addr = mmap((caddr_t)0, maps[i].len,
			PROT_READ, MAP_SHARED, fd, (off_t)0);
		(void)close(fd);
		#endif
		if (maps[i].addr == (caddr_t)-1)
		{
			(void)perror("Memory Map");
			mxFree(buf);
			free_maps(maps,i+1);
			mexErrMsgTxt("Cant map image file.");
		}
		mxFree(buf);
	}


	tmp=mxGetField(ptr,i,"pinfo");
	if (tmp != (mxArray *)0)
	{
		if (mxGetM(tmp) != 3 || (mxGetN(tmp) != 1 && mxGetN(tmp) != maps[i].dim[2]))
		{
			free_maps(maps,i+1);
			mexErrMsgTxt("Wrong sized pinfo.");
		}
		pr = mxGetPr(tmp);
		if (mxGetN(tmp) == 1)
		{
			off = (int)fabs(pr[2]);
			if (off+maps[i].dim[0]*maps[i].dim[1]*maps[i].dim[2]*dsize/8 > maps[i].len)
			{
				free_maps(maps,i+1);
				mexErrMsgTxt("File too small.");
			}
			for(j=0; j<maps[i].dim[2]; j++)
			{
				maps[i].scale[j]  = pr[0];
				maps[i].offset[j] = pr[1];
				maps[i].data[j]   = maps[i].addr+off+j*maps[i].dim[0]*maps[i].dim[1]*dsize/8;
			}
		}
		else
		{
			for(j=0; j<maps[i].dim[2]; j++)
			{
				maps[i].scale[j]  = pr[0+j*3];
				maps[i].offset[j] = pr[1+j*3];
				off = (int)fabs(pr[2+j*3]);
				maps[i].data[j]   = maps[i].addr+off;
				if (off+maps[i].dim[0]*maps[i].dim[1]*dsize/8 > maps[i].len)
				{
					free_maps(maps,i+1);
					mexErrMsgTxt("File too small.");
				}
			}
		}
	}
	else
	{
		if (maps[i].dim[0]*maps[i].dim[1]*maps[i].dim[2]*dsize/8 > maps[i].len)
		{
			free_maps(maps,i+1);
			mexErrMsgTxt("File too small.");
		}
		for(j=0; j<maps[i].dim[2]; j++)
		{
			maps[i].scale[j]  = 1.0;
			maps[i].offset[j] = 0.0;
			maps[i].data[j]   = maps[i].addr+j*maps[i].dim[0]*maps[i].dim[1]*dsize/8;
		}
	}


	tmp=mxGetField(ptr,i,"mat");
	if (tmp != (mxArray *)0)
	{
		if (mxGetM(tmp) != 4 || mxGetN(tmp) != 4)
		{
			free_maps(maps,i+1);
			mexErrMsgTxt("Wrong sized mat.");
		}
		pr = mxGetPr(tmp);
		for(j=0; j<16; j++)
			maps[i].mat[j] = pr[j];
	}
	else
	{
		for(j=0; j<16; j++)
			maps[i].mat[j] = 0.0;
		for(j=0; j<4; j++)
			maps[i].mat[j + j*4] = 1.0;
	}
}

/**************************************************************************/

static MAPTYPE *get_maps_struct(const mxArray *ptr, int *n)
{
	MAPTYPE *maps;
	int num_dims, i;
	const int *dims;

	if (!mxIsStruct(ptr))
	{
		mexErrMsgTxt("Not a structure.");
		return((MAPTYPE *)0);
	}
	num_dims  = mxGetNumberOfDimensions(ptr);
	dims      = mxGetDimensions(ptr);
	*n = 1;
	for(i=0; i<num_dims; i++)
		*n *= dims[i];

	maps = (MAPTYPE *)mxCalloc(*n, sizeof(MAPTYPE));

	if (*n > 0)
	{
		if (mxGetField(ptr,0,"dat") == (mxArray *)0)
			for (i=0; i< *n; i++)
				get_map_file(i, ptr, maps);
		else
			for (i=0; i< *n; i++)
				get_map_dat(i, ptr, maps);
	}
	return(maps);
}

/**************************************************************************/

static MAPTYPE *get_maps_oldstyle(const mxArray *matrix_ptr, int *pn)
{
	MAP *current_map;
	MAPTYPE *maps;
	int n, j;

	if ((!mxIsNumeric(matrix_ptr) || mxIsComplex(matrix_ptr) ||
		mxIsSparse(matrix_ptr) || !mxIsDouble(matrix_ptr) ||
		mxGetM(matrix_ptr) != (sizeof(MAP)+sizeof(double)-1)/sizeof(double)))
	{
		mexErrMsgTxt("Bad image handle dimensions.");
	}
	n = mxGetN(matrix_ptr);

	for (j=0; j<n; j++)
	{
		int xdim, ydim, zdim, datasize;
		current_map = (MAP *)(mxGetPr(matrix_ptr)+
			+ j*((sizeof(MAP)+sizeof(double)-1)/sizeof(double)));
		if (current_map->magic != MAGIC)
		{
			mexErrMsgTxt("Bad magic number in image handle.");
		}
		if (current_map->pid != getpid())
		{
			mexErrMsgTxt("Invalid image handle (from old session).");
		}
		datasize = get_datasize(current_map->dtype);	
		if (datasize == 0)
			mexErrMsgTxt("Bad datatype in image handle.");

		xdim = abs((int)current_map->dim[0]);
		ydim = abs((int)current_map->dim[1]);
		zdim = abs((int)current_map->dim[2]);
		if ((xdim*ydim*zdim*datasize+7)/8 != current_map->len)
		{
			mexErrMsgTxt("Who screwed up the handle??");
		}
	}

	maps = (MAPTYPE *)mxCalloc(n, sizeof(MAPTYPE));
	for (j=0; j<n; j++)
	{
		int jj, t;
		current_map = (MAP *)(mxGetPr(matrix_ptr)
			+ j*((sizeof(MAP)+sizeof(double)-1)/sizeof(double)));

		maps[j].dtype     = current_map->dtype;

		maps[j].dim[0]    = current_map->dim[0];
		maps[j].dim[1]    = current_map->dim[1];
		maps[j].dim[2]    = current_map->dim[2];

		maps[j].data   = (void  **)mxCalloc(maps[j].dim[2],sizeof(void *));
		maps[j].scale  = (double *)mxCalloc(maps[j].dim[2],sizeof(double));
		maps[j].offset = (double *)mxCalloc(maps[j].dim[2],sizeof(double));

		t = maps[j].dim[0]*maps[j].dim[1]*get_datasize(maps[j].dtype)/8;
		for(jj=0; jj<maps[j].dim[2]; jj++)
		{
			maps[j].scale[jj]  = current_map->scale;
			maps[j].offset[jj] = 0.0;
			maps[j].data[jj]   = current_map->data + jj*t;
		}

		for(jj=0; jj<16; jj++)
			maps[j].mat[jj] = 0.0;
		maps[j].mat[0 + 0*4] = current_map->pixdim[0];
		maps[j].mat[1 + 1*4] = current_map->pixdim[1];
		maps[j].mat[2 + 2*4] = current_map->pixdim[2];
		maps[j].mat[3 + 3*4] = 1.0;

		maps[j].addr      = 0;
		maps[j].len       = 0;
	}
	*pn = n;
	return(maps);
}

/**************************************************************************/

static MAPTYPE *get_maps_3dvol(const mxArray *ptr, int *n)
{
	int num_dims, jj, t, dtype;
	const int *dims;
	MAPTYPE *maps;
	unsigned char *dptr;

	if      (mxIsDouble(ptr)) dtype = SPM_DOUBLE;
	else if (mxIsSingle(ptr)) dtype = SPM_FLOAT;
	else if (mxIsInt32 (ptr)) dtype = SPM_SIGNED_INT;
	else if (mxIsUint32(ptr)) dtype = SPM_UNSIGNED_INT;
	else if (mxIsInt16 (ptr)) dtype = SPM_SIGNED_SHORT;
	else if (mxIsUint16(ptr)) dtype = SPM_UNSIGNED_SHORT;
	else if (mxIsInt8  (ptr)) dtype = SPM_SIGNED_CHAR;
	else if (mxIsUint8 (ptr)) dtype = SPM_UNSIGNED_CHAR;
	else mexErrMsgTxt("Unknown volume datatype.");

	maps = (MAPTYPE *)mxCalloc(1, sizeof(MAPTYPE));

	num_dims = mxGetNumberOfDimensions(ptr);
	if (num_dims > 3)
	{
		mexErrMsgTxt("Too many dimensions.");
	}

	dims     = mxGetDimensions(ptr);
	for(jj=0; jj<num_dims; jj++)
		maps->dim[jj]=dims[jj];
	for(jj=num_dims; jj<3; jj++)
		maps->dim[jj]=1;

	for(jj=0; jj<16; jj++)
		maps->mat[jj] = 0.0;
	for(jj=0; jj<4; jj++)
		maps->mat[jj + jj*4] = 1.0;

	maps->dtype  = dtype;

	maps->data   = (void  **)mxCalloc(maps->dim[2],sizeof(void *));
	maps->scale  = (double *)mxCalloc(maps->dim[2],sizeof(double));
	maps->offset = (double *)mxCalloc(maps->dim[2],sizeof(double));

	t     = maps->dim[0]*maps->dim[1]*get_datasize(maps->dtype)/8;
	dptr   = (unsigned char *)mxGetPr(ptr);

	for(jj=0; jj<maps->dim[2]; jj++)
	{
		maps->scale[jj]  = 1.0;
		maps->offset[jj] = 0.0;
		maps->data[jj]   = &(dptr[jj*t]);
	}

	maps->addr      = 0;
	maps->len       = 0;

	*n = 1;
	return(maps);
}

/**************************************************************************/

MAPTYPE *get_maps(const mxArray *ptr, int *n)
{
	if (mxIsStruct(ptr))
		return(get_maps_struct(ptr, n));
	/* else if (mxIsNumeric(ptr) && !mxIsComplex(ptr) &&
		!mxIsSparse(ptr) && mxIsDouble(ptr) &&
		mxGetM(ptr) == (sizeof(MAP)+sizeof(double)-1)/sizeof(double))
		return(get_maps_oldstyle(ptr, n)); */
	else if (mxGetNumberOfDimensions(ptr) <= 3 &&
		mxIsNumeric(ptr) && !mxIsComplex(ptr) &&
		!mxIsSparse(ptr))
		return(get_maps_3dvol(ptr, n));
	else
		mexErrMsgTxt("What do I do with this?");
	return((MAPTYPE *)0);
}


/**************************************************************************/

void voxdim(MAPTYPE *map, double vdim[3])
{
	int i, j;
	double t;
	for(j=0; j<3; j++)
	{
		t=0.0;
		for(i=0; i<3; i++)
			t += map->mat[i+j*4]*map->mat[i+j*4];
		vdim[j] = sqrt(t);
	}
}

/**************************************************************************/

int get_dtype(const mxArray *ptr)
{
	mxArray *tmp;
	double *pr;

	if (!mxIsStruct(ptr))
	{
		mexErrMsgTxt("Not a structure.");
		return(0);
	}
	tmp=mxGetField(ptr,0,"dim");
	if (tmp == (mxArray *)0)
	{
		mexErrMsgTxt("Cant find dim.");
	}
	if (mxGetM(tmp)*mxGetN(tmp) != 4)
	{
		mexErrMsgTxt("Wrong sized dim.");
	}
	pr = mxGetPr(tmp);

	return((int)fabs(pr[3]));
}

/**************************************************************************/
