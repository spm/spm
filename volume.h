#ifndef volume_hdr
#define volume_hdr

#ifndef lint
static char hdr_sccsid[]="%W% John Ashburner %E%";
#endif
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include "cmex.h"

#define MAGIC 110494
typedef struct maptype
{
	double xdim, ydim, zdim;
	double xpixdim, ypixdim, zpixdim;
	double scalefactor;
	int datatype;
	int off;
	caddr_t map;
	unsigned char *data;
	size_t len;
	int magic;
	int prot;
	int flags;
/*
	int fd;
*/
	pid_t pid;
}	MAPTYPE;

/* a few datatypes nicked from ANALYZE */
#define UNSIGNED_CHAR	2
#define SIGNED_SHORT	4
#define SIGNED_INT	8
#define FLOAT	16
#define DOUBLE	64

int get_datasize(type)
int type;
{
	if (type == UNSIGNED_CHAR) return(8);
	if (type == SIGNED_SHORT) return(16);
	if (type == SIGNED_INT) return(32);
	if (type == FLOAT) return(32);
	if (type == DOUBLE) return(64);
	return(0);
}

MAPTYPE *get_map(matrix_ptr)
Matrix *matrix_ptr;
{
	MAPTYPE *map;
	static struct stat stbuf;
	int xdim, ydim, zdim, datasize;

	if ((!mxIsNumeric(matrix_ptr) || mxIsComplex(matrix_ptr) ||
		!mxIsFull(matrix_ptr) || !mxIsDouble(matrix_ptr) ||
		mxGetN(matrix_ptr) != 1 ||
		mxGetM(matrix_ptr) != (sizeof(MAPTYPE)+sizeof(double)-1)/sizeof(double)))
	{
		mexErrMsgTxt("Bad image handle dimensions.");
	}
	map = (MAPTYPE *)mxGetPr(matrix_ptr);
	if (map->magic != MAGIC)
	{
		mexErrMsgTxt("Bad magic number in image handle.");
	}
	if (map->pid != getpid())
	{
		mexErrMsgTxt("Invalid image handle (from old session).");
	}
	datasize = get_datasize(map->datatype);	
	if (datasize == 0)
		mexErrMsgTxt("Bad datatype in image handle.");

	xdim = abs(nint(map->xdim));
	ydim = abs(nint(map->ydim));
	zdim = abs(nint(map->zdim));
	if ((xdim*ydim*zdim*datasize+7)/8 != map->len)
	{
		mexErrMsgTxt("Who screwed up the handle??");
	}
/*
	if (fstat(map->fd, &stbuf) == -1)
	{
		mexErrMsgTxt("Cant stat image file.");
	}
	if (stbuf.st_size < map->off+map->len)
	{
		mexErrMsgTxt("Image file too small.");
	}
*/
	return(map);
}
#endif
