#ifndef lint
static char sccsid[]="%W% John Ashburner %E%";
#endif

#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include "cmex.h"
#include "volume.h"

#ifdef __STDC__
void mexFunction(int nlhs, Matrix *plhs[], int nrhs, Matrix *prhs[])
#else
mexFunction(nlhs, plhs, nrhs, prhs)
int nlhs, nrhs;
Matrix *plhs[], *prhs[];
#endif
{
	char *str;
	int k,stlen, datasize;
	MAPTYPE *map;
	double *ptr;
	static struct stat stbuf;

	if (nrhs != 2 || nlhs > 1)
		mexErrMsgTxt("Inappropriate usage.");

	/* get filename */
	if (!mxIsString(prhs[0]))
		mexErrMsgTxt("filename should be a string");
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
		!mxIsFull(prhs[1]) || !mxIsDouble(prhs[1]))
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

	plhs[0] = mxCreateFull((sizeof(MAPTYPE)+sizeof(double)-1)/sizeof(double),1,REAL);
	map = (MAPTYPE *)mxGetPr(plhs[0]);

	map->magic = MAGIC;
	map->xdim = abs(nint(ptr[0]));
	map->ydim = abs(nint(ptr[1]));
	map->zdim = abs(nint(ptr[2]));
	map->xpixdim = ptr[3];
	map->ypixdim = ptr[4];
	map->zpixdim = ptr[5];
	map->scalefactor = ptr[6];
	map->datatype = abs(nint(ptr[7]));
	map->off = abs(nint(ptr[8]));
	map->pid = getpid();

	if (map->datatype == UNSIGNED_CHAR)
	{
		datasize = 8;
	}
	else if (map->datatype == SIGNED_SHORT)
	{
		datasize = 16;
	}
	else
	{
		mxFree(str);
		mexErrMsgTxt("Unrecognised datatype.");
	}
	map->len = ((int)(map->xdim*map->ydim*map->zdim)*datasize+7)/8;
	map->prot = PROT_READ;
	map->flags = MAP_SHARED;

	if ((map->fd = open(str, O_RDONLY)) == -1)
	{
		mxFree(str);
		mexErrMsgTxt("Cant open image file.");
	}

	if (fstat(map->fd, &stbuf) == -1)
	{
		close(map->fd);
		mxFree(str);
		mexErrMsgTxt("Cant stat image file.");
	}
	if (stbuf.st_size < map->off+map->len)
	{
		close(map->fd);
		mxFree(str);
		mexErrMsgTxt("Image file too small.");
	}

	map->map = mmap((caddr_t)0, map->len, map->prot, map->flags, map->fd, map->off);
	if (map->map == (caddr_t)-1)
	{
		close(map->fd);
		mxFree(str);
		mexErrMsgTxt("Cant map image file.");
	}
	mxFree(str);
}
