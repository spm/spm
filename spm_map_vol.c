#ifndef lint
static char sccsid[]="%W% John Ashburner %E%";
#endif

#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include "mex.h"
#include "spm_map.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	char *str, errstr[2048];
	int k,stlen, datasize, fd;
	MAP *map;
	double *ptr;
	static struct stat stbuf;

	if (nrhs != 2 || nlhs > 1)
		mexErrMsgTxt("Inappropriate usage.");

	/* get filename 
	if (!mxIsString(prhs[0]))
		mexErrMsgTxt("filename should be a string"); */
	stlen = mxGetN(prhs[0]);
	str = (char *)mxCalloc(stlen+1, sizeof(char));
	mxGetString(prhs[0],str,stlen+1);

	/* delete white space */
	for(k=0; k<stlen; k++)
		if (str[k] == ' ')
		{
			str[k] = '\0';
			break;
		}

	/* get options */
	if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) ||
		mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1]))
	{
		mxFree(str);
		mexErrMsgTxt("Options vector must be numeric, real, full and double.");
	}
	if (mxGetM(prhs[1])* mxGetN(prhs[1]) != 9)
	{
		mxFree(str);
		mexErrMsgTxt("Options vector is wrong size.");
	}
	ptr = mxGetPr(prhs[1]);

	plhs[0] = mxCreateDoubleMatrix((sizeof(MAP)+sizeof(double)-1)/sizeof(double),1,mxREAL);
	map = (MAP *)mxGetPr(plhs[0]);

	map->magic = MAGIC;
	map->dim[0] = abs((int)ptr[0]);
	map->dim[1] = abs((int)ptr[1]);
	map->dim[2] = abs((int)ptr[2]);
	map->pixdim[0] = ptr[3];
	map->pixdim[1] = ptr[4];
	map->pixdim[2] = ptr[5];
	map->scale = ptr[6];
	map->dtype = abs((int)ptr[7]);
	map->off = abs((int)ptr[8]);
	map->pid = getpid();
	if (map->dim[0] < 1 || map->dim[1] < 1 || map->dim[2] < 1)
		mexErrMsgTxt("Dimensions too small.");

	datasize = get_datasize(map->dtype);
	if (datasize == 0)
	{
		mxFree(str);
		mexErrMsgTxt("Unrecognised datatype.");
	}
	map->len = ((int)(map->dim[0]*map->dim[1]*map->dim[2])*datasize+7)/8;
	map->prot = PROT_READ;
	map->flags = MAP_SHARED;

	if ((fd = open(str, O_RDONLY)) == -1)
	{
		(void)sprintf(errstr,"Cant open image file (%s).", str);
		mxFree(str);
		mexErrMsgTxt(errstr);
	}

	if (fstat(fd, &stbuf) == -1)
	{
		(void)close(fd);
		(void)sprintf(errstr,"Cant stat image file (%s).", str);
		mxFree(str);
		mexErrMsgTxt(errstr);
	}
	if (stbuf.st_size < map->off+map->len)
	{
		(void)close(fd);
		(void)sprintf(errstr,"Image file too small (%s).", str);
		mxFree(str);
		mexErrMsgTxt(errstr);
	}

	map->map = mmap((caddr_t)0, map->len+map->off, map->prot, map->flags, fd, (off_t)0);
	if (map->map == (caddr_t)-1)
	{
		(void)perror("Memory Map");
		(void)close(fd);
		(void)sprintf(errstr,"Cant map image file (%s).", str);
		mxFree(str);
		mexErrMsgTxt(errstr);
	}
	map->data = (unsigned char *)(map->map) + map->off;
	mxFree(str);
	(void)close(fd);
}
