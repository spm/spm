#ifndef lint
static char sccsid[]="%W% John Ashburner %E%";
#endif

#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include "cmex.h"
#include "dbh.h"
#define MAGIC 170372
typedef struct mayomaptype
{
	int magic;
	struct dsr hdr;
	caddr_t map;
	size_t len;
	int prot;
	int flags;
}	MAYOMAPTYPE;


#ifdef __STDC__
void mexFunction(int nlhs, Matrix *plhs[], int nrhs, Matrix *prhs[])
#else
mexFunction(nlhs, plhs, nrhs, prhs)
int nlhs, nrhs;
Matrix *plhs[], *prhs[];
#endif
{
	char *str;
	int k,stlen, fd;
	MAYOMAPTYPE *mayomap;
	struct image_dimension *dime;

	if (nrhs != 1 || nlhs > 1)
		mexErrMsgTxt("inappropriate usage");

	/* get filename */
	if (!mxIsString(prhs[0]))
		mexErrMsgTxt("filename should be a string");
	stlen = mxGetN(prhs[0]);
	str = (char *)mxCalloc(stlen+1+4, sizeof(char));
	mxGetString(prhs[0],str,stlen+1);

	/* delete white space */
	for(k=0; k<stlen; k++)
		if (str[k] == ' ')
		{
			str[k] = '\0';
			stlen = k;
			break;
		}

	/* check for suffixes */
	if (!strcmp(str+stlen-4,".img") || !strcmp(str+stlen-4,".hdr"))
		stlen -= 4;

	/* load .hdr file */
	sprintf(str+stlen, "%s", ".hdr");
	if ((fd = open(str, O_RDONLY)) == 0)
	{
		mxFree(str);
		mexErrMsgTxt("cant open .hdr file");
	}
	plhs[0] = mxCreateFull((sizeof(MAYOMAPTYPE)+sizeof(double)-1)/sizeof(double),1,REAL);
	mayomap = (MAYOMAPTYPE *)mxGetPr(plhs[0]);
	mayomap->magic = MAGIC;
	if (read(fd, (char *)(&(mayomap->hdr)), sizeof(struct dsr)) != sizeof(struct dsr))
	{
		close(fd);
		mxFree(str);
		mexErrMsgTxt("cant read .hdr file");
	}
	close(fd);

	/* map image file */
	sprintf(str+stlen, "%s", ".img");
	if ((fd = open(str, O_RDONLY)) == 0)
	{
		mexErrMsgTxt("cant open .img file");
	}
	dime = &(mayomap->hdr.dime);
	mayomap->len = (dime->dim[1]*dime->dim[2]*dime->dim[3]*dime->bitpix+7)/8;
	mayomap->prot = PROT_READ;
	mayomap->flags = MAP_SHARED;

	mayomap->map = mmap((caddr_t)0, mayomap->len,
		mayomap->prot , mayomap->flags, fd, (off_t)0);
	if (mayomap->map == (caddr_t)-1)
	{
		close(fd);
		mxFree(str);
		mexErrMsgTxt("cant map .img file");
	}
	close(fd);

	mxFree(str);
}
