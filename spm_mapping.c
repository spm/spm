#ifndef lint
static char sccsid[]="%W% John Ashburner %E%";
#endif

#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <math.h>
#include "mex.h"
#include "spm_vol_utils.h"
#include "spm_map.h"


int get_datasize(int type)
{
	if (type == UNSIGNED_CHAR || type == SIGNED_CHAR) return(8);
	if (type == SIGNED_SHORT || type == SIGNED_SHORT_S) return(16);
	if (type == UNSIGNED_SHORT || type == UNSIGNED_SHORT_S) return(16);
	if (type == SIGNED_INT || type == SIGNED_INT_S) return(32);
	if (type == UNSIGNED_INT || type == UNSIGNED_INT_S) return(32);
	if (type == FLOAT || type == FLOAT_S) return(32);
	if (type == DOUBLE || type == DOUBLE_S) return(64);
	return(0);
}

void free_maps(MAPTYPE *maps, int n)
{
	int j;
	for(j=0; j<n; j++)
	{
		if (maps[j].addr)
		{
			(void)munmap(maps[j].addr, maps[j].len);
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

MAPTYPE *get_maps_struct(const mxArray *ptr, int *n)
{
	MAPTYPE *maps;
	int num_dims, i,j;
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


	for (i=0; i<*n; i++)
	{
		mxArray *tmp;
		double *pr;
		int dsize, off;

		tmp=mxGetField(ptr,i,"dim");
		if (tmp == (mxArray *)0)
		{
			free_maps(maps,i+1);
			mexErrMsgTxt("Cant find dim.");
		}
		if (mxGetM(tmp)*mxGetN(tmp) != 4)
		{
			free_maps(maps,i+1);
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
			maps[i].addr = mmap((caddr_t)0, maps[i].len,
				PROT_READ, MAP_SHARED, fd, (off_t)0);
			if (maps[i].addr == (caddr_t)-1)
			{
				(void)perror("Memory Map");
				(void)close(fd);
				mxFree(buf);
				free_maps(maps,i+1);
				mexErrMsgTxt("Cant map image file.");
			}
			(void)close(fd);
			mxFree(buf);
		}


		tmp=mxGetField(ptr,i,"mat");
		if (tmp == (mxArray *)0)
		{
			free_maps(maps,i+1);
			mexErrMsgTxt("Cant find mat.");
		}
		if (mxGetM(tmp) != 4 || mxGetN(tmp) != 4)
		{
			free_maps(maps,i+1);
			mexErrMsgTxt("Wrong sized mat.");
		}
		pr = mxGetPr(tmp);
		for(j=0; j<16; j++)
		{
			maps[i].mat[j] = pr[j];
		}

		tmp=mxGetField(ptr,i,"pinfo");
		if (tmp == (mxArray *)0)
		{
			free_maps(maps,i+1);
			mexErrMsgTxt("Cant find pinfo.");
		}
		if (mxGetM(tmp) != 3 || (mxGetN(tmp) != 1 && mxGetN(tmp) != maps[i].dim[2]))
		{
			free_maps(maps,i+1);
			mexErrMsgTxt("Wrong sized pinfo.");
		}
		pr = mxGetPr(tmp);

		if (mxGetN(tmp) == 1)
		{
			int j;
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
			int j;
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
	return(maps);
}

MAPTYPE *get_maps_oldstyle(const mxArray *matrix_ptr, int *pn)
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

MAPTYPE *get_maps_3dvol(const mxArray *ptr, int *n)
{
	int num_dims, jj, t;
	const int *dims;
	MAPTYPE *maps;

	maps = (MAPTYPE *)mxCalloc(1, sizeof(MAPTYPE));

	num_dims = mxGetNumberOfDimensions(ptr);
	dims     = mxGetDimensions(ptr);
	for(jj=0; jj<num_dims; jj++)
		maps->dim[jj]=dims[jj];
	for(jj=num_dims; jj<3; jj++)
		maps->dim[jj]=1;

	for(jj=0; jj<16; jj++)
		maps->mat[jj] = 0.0;
	for(jj=0; jj<4; jj++)
		maps->mat[jj + jj*4] = 1.0;

	maps->dtype     = DOUBLE;

	maps->data   = (void  **)mxCalloc(maps->dim[2],sizeof(void *));
	maps->scale  = (double *)mxCalloc(maps->dim[2],sizeof(double));
	maps->offset = (double *)mxCalloc(maps->dim[2],sizeof(double));

	t = maps->dim[0]*maps->dim[1];
	for(jj=0; jj<maps->dim[2]; jj++)
	{
		maps->scale[jj]  = 1.0;
		maps->offset[jj] = 0.0;
		maps->data[jj]   = &(mxGetPr(ptr)[jj*t]);
	}

	maps->addr      = 0;
	maps->len       = 0;

	*n = 1;
	return(maps);
}

MAPTYPE *get_maps(const mxArray *ptr, int *n)
{
	if (mxIsStruct(ptr))
		return(get_maps_struct(ptr, n));
	else if (mxIsNumeric(ptr) && !mxIsComplex(ptr) &&
		!mxIsSparse(ptr) && mxIsDouble(ptr) &&
		mxGetM(ptr) == (sizeof(MAP)+sizeof(double)-1)/sizeof(double))
		return(get_maps_oldstyle(ptr, n));
	else if (mxGetNumberOfDimensions(ptr) <= 3 &&
		mxIsNumeric(ptr) && !mxIsComplex(ptr) &&
		!mxIsSparse(ptr) && mxIsDouble(ptr))
		return(get_maps_3dvol(ptr, n));
	else
		mexErrMsgTxt("What do I do with this?");
	return((MAPTYPE *)0);
}


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
