#ifndef lint
static char sccsid[]="%W% (c) John Ashburner %E%";
#endif
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/param.h>
#include <sys/stat.h>
#include <dirent.h>

#include "mex.h"


slowsort(argc, argv)
char **argv;
int argc;
{
	int i, j;
	for(i=argc; i>1; i--)
		for(j=1; j<i; j++)
		{
			if(strcmp(argv[j], argv[j-1]) < 0)
			{
				char *tempp;
				tempp = argv[j];
				argv[j] = argv[j-1];
				argv[j-1] = tempp;
			}
		}
}

wildcard(attempt, actual)
char *attempt;
char *actual;
{
	char *p1, *p2;

	for(p1=attempt, p2=actual;
		((*p1 == *p2) || (*p1 == '?')) && ((*p1 != '\0') && (*p2 != '\0'));
		p1++, p2++)
		;
	if ((*p1 == *p2) || ((*p1 == '?') && (*p2 != '\0')))
		return 1;
	if (*p1 != '*')
		return 0;

	for (;*p1 == '*';)
		p1++;
	for (;(wildcard(p1, p2) == 0) && (*p2 != '\0');)
		p2++;
	if ((*p1 == *p2) || ((*p1 == '?') && (*p2 != '\0')))
		return 1;
	return 0;
}

void list2mat(m,n,list,ptr)
int m,n;
char *list[];
Matrix **ptr;
{
	char *buf2;
	int i;
	buf2 = (char *)mxCalloc(m*n+1,1);
	for (i=0; i<m; i++)
	{
		int j = 0;
		while(list[i][j])
		{
			buf2[i+j*m]=list[i][j];
			j++;
		}
		while (j < n)
		{
			buf2[i+j*m]=' ';
			j++;
		}
	}
	buf2[m*n] = 0;
	*ptr = mxCreateString(buf2);
	mxSetM(*ptr,m);
	mxSetN(*ptr,n);
	(void)mxFree(buf2);
}


#ifdef __STDC__
void mexFunction(int nlhs, Matrix *plhs[], int nrhs, Matrix *prhs[])
#else
mexFunction(nlhs, plhs, nrhs, prhs)
int nlhs, nrhs;
Matrix *plhs[], *prhs[];
#endif
{
	int ndirs = 0, nfiles = 0, len, uid, gids[128], ngids, maxdlen = 0, maxflen = 0;
	DIR *dirp;
	struct dirent *dp;
	char *filenames[1024], *directories[1024], *ptr, *buf = (char *)0, *bufp, *fullpathname, *filter;
	static struct stat stbuf;

	if ((nrhs != 2) || (nlhs != 2))
		mexErrMsgTxt("Incorrect Usage.");

	if (!mxIsString(prhs[0]) || !mxIsString(prhs[1]))
		mexErrMsgTxt("Arguments must be strings.");

	if (mxGetM(prhs[0]) != 1 || mxGetM(prhs[1]) != 1)
		mexErrMsgTxt("Arguments must have one row.");

	len = mxGetN(prhs[0]);
	fullpathname = (char *)mxCalloc(len+1+MAXNAMLEN, sizeof(char));
	mxGetString(prhs[0],fullpathname,len+1);

	len = mxGetN(prhs[1]);
	filter = (char *)mxCalloc(len+1, sizeof(char));
	mxGetString(prhs[1],filter,len+1);

	uid = getuid();
	ngids = getgroups(128,gids);

	if ((stat(fullpathname, &stbuf) != -1) && (dirp = opendir(fullpathname)))
	{
		buf = (char *)mxCalloc((int)stbuf.st_size, 1);
		bufp = buf;
		ptr = fullpathname + strlen(fullpathname);
		(void)strcat(fullpathname, "/");
		ptr++;

		for (dp = readdir(dirp); dp != NULL; dp = readdir(dirp))
		{
			*ptr = 0;
			(void)strcpy(ptr, dp->d_name);
			if (stat(fullpathname, &stbuf) != -1)
			{
				int g;
				int mask = 0007;
				if ((stbuf.st_mode & S_IFMT) == S_IFDIR
					|| ((stbuf.st_mode & S_IFMT) == S_IFREG && wildcard(filter,dp->d_name)))
				{
					if (stbuf.st_uid == uid)
						mask = 0700; /* user */
					else
						for(g=0; g<ngids; g++)
							if (gids[g] == stbuf.st_gid)
								mask = 0070; /* group */

					if ((stbuf.st_mode & S_IFMT) == S_IFDIR && (mask & 0555 & stbuf.st_mode))
					{
						if (ndirs == 1024)
							mexErrMsgTxt("Too many subdirectories.");

						len = strlen(dp->d_name);
						if (len > maxdlen) maxdlen = len;
						directories[ndirs] = bufp;
						bufp += (len+1);
						(void)strcpy(directories[ndirs], dp->d_name);
						ndirs++;
					}
					else if ((stbuf.st_mode & S_IFMT) == S_IFREG &&  (mask & 0444 & stbuf.st_mode))
					{
						if (nfiles == 1024)
							mexErrMsgTxt("Too many files match.");
						len = strlen(dp->d_name);
						if (len > maxflen) maxflen = len;

						filenames[nfiles] = bufp;
						bufp += (len+1);
						(void)strcpy(filenames[nfiles], dp->d_name);
						nfiles++;
					}
				}
			}
		}
		(void)closedir(dirp);
	}
	slowsort(nfiles, filenames);
	list2mat(nfiles,maxflen,filenames,&plhs[0]);

	slowsort(ndirs, directories);
	list2mat(ndirs,maxdlen,directories,&plhs[1]);

	if (buf) (void)mxFree(buf);
}

