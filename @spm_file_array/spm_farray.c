/*
 * $Id: spm_farray.c 253 2005-10-13 15:31:34Z guillaume $
 */

#include <math.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <sys/types.h>

#ifdef SPM_WIN32
#include <windows.h>
#include <memory.h>
typedef char *caddr_t;
#else
#include <unistd.h>
#endif

#include "mex.h"
#include "spm_mapping.h"
#include "spm_datatypes.h"

#define MXDIMS 256

int icumprod[MXDIMS], ocumprod[MXDIMS];

typedef struct mtype {
	int ndim;
	int dim[MXDIMS];
	int dtype;
	caddr_t addr;
	size_t len;
	void *data;
} MTYPE;

void do_unmap_file(MTYPE *map)
{
	if (map->addr)
	{
#ifdef SPM_WIN32
		(void)unmap_file(map->addr);
#else
		(void)munmap(map->addr, map->len);
#endif
		map->addr=0;
	}
}

void do_map_file(const mxArray *ptr, MTYPE *map)
{
	int i;
	double *pr;
	mxArray *tmp;
	int siz;
	if (!mxIsStruct(ptr))
	{
		mexErrMsgTxt("Not a structure.");
	}

	tmp=mxGetField(ptr,0,"dtype");
	if (tmp == (mxArray *)0)
	{
		mexErrMsgTxt("Cant find dtype.");
	}
	if (mxGetM(tmp)*mxGetN(tmp) != 1)
	{
		mexErrMsgTxt("Datatype should be a scaler.");
	}
	map->dtype = mxGetPr(tmp)[0];
        if (map->dtype != SPM_UNSIGNED_CHAR && map->dtype != SPM_SIGNED_CHAR && 
        	map->dtype != SPM_SIGNED_SHORT && map->dtype != SPM_UNSIGNED_SHORT &&
        	map->dtype != SPM_SIGNED_INT && map->dtype != SPM_UNSIGNED_INT &&
		map->dtype != SPM_FLOAT && map->dtype != SPM_DOUBLE)
	{
		mexErrMsgTxt("Bad datatype.");
	}

	tmp      = mxGetField(ptr,0,"dim");
	if (tmp == (mxArray *)0)
	{
		mexErrMsgTxt("Cant find dim.");
	}
	map->ndim = mxGetM(tmp)*mxGetN(tmp);
	if (map->ndim >= MXDIMS)
	{
		mexErrMsgTxt("Too many dimensions.");
	}
	pr       = mxGetPr(tmp);
	siz      = 1;
	for(i=0; i<map->ndim; i++)
	{
		map->dim[i] = (int)fabs(pr[i]);
		siz = siz*map->dim[i];
	}

	tmp=mxGetField(ptr,0,"fname");
	if (tmp == (mxArray *)0)
	{
		mexErrMsgTxt("Cant find fname.");
	}
	if (mxIsChar(tmp))
	{
		int buflen;
		char *buf;
		int fd;
		struct stat stbuf;
		buflen = mxGetN(tmp)*mxGetM(tmp)+1;
		buf    = mxCalloc(buflen,sizeof(char));
		if (mxGetString(tmp,buf,buflen))
		{
			mxFree(buf);
			mexErrMsgTxt("Cant get filename.");
		}
		if ((fd = open(buf, O_RDONLY)) == -1)
		{
			mxFree(buf);
			mexErrMsgTxt("Cant open file.");
		}
		if (fstat(fd, &stbuf) == -1)
		{
			(void)close(fd);
			mxFree(buf);
			mexErrMsgTxt("Cant stat image file.");
		}
		map->len = stbuf.st_size;
		if (map->len < siz)
		{
			(void)close(fd);
			mxFree(buf);
			mexErrMsgTxt("File too small.");
		}
#ifdef SPM_WIN32
		(void)close(fd);
		map->addr = map_file(buf, (caddr_t)0, map->len,
			PROT_READ, MAP_SHARED, (off_t)0);
#else
		map->addr = mmap((caddr_t)0, map->len,
			PROT_READ, MAP_SHARED, fd, (off_t)0);
		(void)close(fd);
#endif
		if (map->addr == (caddr_t)-1)
		{
			(void)perror("Memory Map");
			mxFree(buf);
			mexErrMsgTxt("Cant map image file.");
		}
		mxFree(buf);
	}
	map->data = (void *)map->addr;
}

void get_values_uint8(int ndim, int idim[], unsigned char idat[], int odim[], double odat[], int *iptr[])
{
	int i;
	if (ndim == 0)
	{
		for(i=0; i<odim[0]; i++)
			odat[i] = (double)idat[iptr[0][i]-1];
	}
	else
	{
		for(i=0; i<odim[ndim]; i++)
			get_values_uint8(ndim-1, idim, idat+icumprod[ndim]*(iptr[ndim][i]-1),
				odim, odat+ocumprod[ndim]*i, iptr);
	}
}

void get_values_int8(int ndim, int idim[], char idat[], int odim[], double odat[], int *iptr[])
{
	int i;
	if (ndim == 0)
	{
		for(i=0; i<odim[0]; i++)
			odat[i] = (double)idat[iptr[0][i]-1];
	}
	else
	{
		for(i=0; i<odim[ndim]; i++)
			get_values_int8(ndim-1, idim, idat+icumprod[ndim]*(iptr[ndim][i]-1),
				odim, odat+ocumprod[ndim]*i, iptr);
	}
}

void get_values_uint16(int ndim, int idim[], unsigned short idat[], int odim[], double odat[], int *iptr[])
{
	int i;
	if (ndim == 0)
	{
		for(i=0; i<odim[0]; i++)
			odat[i] = (double)idat[iptr[0][i]-1];
	}
	else
	{
		for(i=0; i<odim[ndim]; i++)
			get_values_uint16(ndim-1, idim, idat+icumprod[ndim]*(iptr[ndim][i]-1),
				odim, odat+ocumprod[ndim]*i, iptr);
	}
}

void get_values_int16(int ndim, int idim[], short idat[], int odim[], double odat[], int *iptr[])
{
	int i;
	if (ndim == 0)
	{
		for(i=0; i<odim[0]; i++)
			odat[i] = (double)idat[iptr[0][i]-1];
	}
	else
	{
		for(i=0; i<odim[ndim]; i++)
			get_values_int16(ndim-1, idim, idat+icumprod[ndim]*(iptr[ndim][i]-1),
				odim, odat+ocumprod[ndim]*i, iptr);
	}
}

void get_values_uint32(int ndim, int idim[], unsigned int idat[], int odim[], double odat[], int *iptr[])
{
	int i;
	if (ndim == 0)
	{
		for(i=0; i<odim[0]; i++)
			odat[i] = (double)idat[iptr[0][i]-1];
	}
	else
	{
		for(i=0; i<odim[ndim]; i++)
			get_values_uint32(ndim-1, idim, idat+icumprod[ndim]*(iptr[ndim][i]-1),
				odim, odat+ocumprod[ndim]*i, iptr);
	}
}

void get_values_int32(int ndim, int idim[], int idat[], int odim[], double odat[], int *iptr[])
{
	int i;
	if (ndim == 0)
	{
		for(i=0; i<odim[0]; i++)
			odat[i] = (double)idat[iptr[0][i]-1];
	}
	else
	{
		for(i=0; i<odim[ndim]; i++)
			get_values_int32(ndim-1, idim, idat+icumprod[ndim]*(iptr[ndim][i]-1),
				odim, odat+ocumprod[ndim]*i, iptr);
	}
}

void get_values_float(int ndim, int idim[], float idat[], int odim[], double odat[], int *iptr[])
{
	int i;
	if (ndim == 0)
	{
		for(i=0; i<odim[0]; i++)
			odat[i] = (double)idat[iptr[0][i]-1];
	}
	else
	{
		for(i=0; i<odim[ndim]; i++)
			get_values_float(ndim-1, idim, idat+icumprod[ndim]*(iptr[ndim][i]-1),
				odim, odat+ocumprod[ndim]*i, iptr);
	}
}

void get_values_double(int ndim, int idim[], double idat[], int odim[], double odat[], int *iptr[])
{
	int i;
	if (ndim == 0)
	{
		for(i=0; i<odim[0]; i++)
			odat[i] = (double)idat[iptr[0][i]-1];
	}
	else
	{
		for(i=0; i<odim[ndim]; i++)
			get_values_double(ndim-1, idim, idat+icumprod[ndim]*(iptr[ndim][i]-1),
				odim, odat+ocumprod[ndim]*i, iptr);
	}
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	MTYPE map;
	double *odat;
	void *idat;
	int i;
	int **iptr, *odim, *idim, ndim;
	int one[1];
	one[0] = 1;

	if (nrhs<2 || nlhs>1) mexErrMsgTxt("Incorrect usage.");

	do_map_file(prhs[0], &map);

	/*
	ndim = mxGetNumberOfDimensions(prhs[0]);
	idim = mxGetDimensions(prhs[0]);
	idat = mxGetData(prhs[0]);
	*/
	ndim = map.ndim;
	idim = map.dim;
	idat = map.data;

	if (ndim >= MXDIMS)
		mexErrMsgTxt("Too many dimensions.");

	if (nrhs > ndim+1) mexErrMsgTxt("Index exceeds matrix dimensions (1).");

	iptr = (int **)mxCalloc(ndim,sizeof(int *));
	odim = (int *)mxCalloc(ndim,sizeof(int));

	for(i=0;i<nrhs-1; i++)
	{
		int j;
		if (!mxIsNumeric(prhs[i+1]) || !mxIsInt32(prhs[i+1]) || mxIsComplex(prhs[i+1]))
			mexErrMsgTxt("Indices must be int32.");
		odim[i] = mxGetM(prhs[i+1])*mxGetN(prhs[i+1]);
		iptr[i] = (int *)mxGetPr(prhs[i+1]);
		for(j=0; j<odim[i]; j++)
			if (iptr[i][j]<1 || iptr[i][j]>idim[i])
				mexErrMsgTxt("Index exceeds matrix dimensions (2).");
	}
	for(i=nrhs-1; i<ndim; i++)
	{
		if (idim[i] == 0)
			mexErrMsgTxt("Index exceeds matrix dimensions (3).");
		odim[i] = 1;
		iptr[i] = one;
	}

	icumprod[0] = 1;
	ocumprod[0] = 1;
	for(i=0; i<ndim; i++)
	{
		icumprod[i+1] = icumprod[i]*idim[i];
		ocumprod[i+1] = ocumprod[i]*odim[i];
	}

	plhs[0] = mxCreateNumericArray(ndim,odim,mxDOUBLE_CLASS,mxREAL);
	odat    = mxGetData(plhs[0]);

	if (map.dtype == SPM_UNSIGNED_CHAR)
		get_values_uint8(ndim-1, idim, (unsigned char *)idat, odim, odat, iptr);
	else if(map.dtype == SPM_SIGNED_CHAR)
		get_values_int8(ndim-1, idim, (char *)idat, odim, odat, iptr);
	else if(map.dtype == SPM_UNSIGNED_SHORT)
		get_values_uint16(ndim-1, idim, (unsigned short *)idat, odim, odat, iptr);
	else if(map.dtype == SPM_SIGNED_SHORT)
		get_values_int16(ndim-1, idim, (short *)idat, odim, odat, iptr);
	else if(map.dtype == SPM_UNSIGNED_INT)
		get_values_uint32(ndim-1, idim, (unsigned int *)idat, odim, odat, iptr);
	else if(map.dtype == SPM_SIGNED_INT)
		get_values_int32(ndim-1, idim, (int *)idat, odim, odat, iptr);
	else if(map.dtype == SPM_FLOAT)
		get_values_float(ndim-1, idim, (float *)idat, odim, odat, iptr);
	else if(map.dtype == SPM_DOUBLE)
		get_values_double(ndim-1, idim, (double *)idat, odim, odat, iptr);
	else
		mexErrMsgTxt("Bad datatype.");


	mxFree((char *)iptr);
	mxFree(odim);
	do_unmap_file(&map);
}
