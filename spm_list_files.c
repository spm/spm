#ifndef lint
static char sccsid[]="%W% (c) John Ashburner %E%";
#endif
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/param.h>
#include <sys/stat.h>
#include <dirent.h>
#include <unistd.h>

#include "mex.h"

#define MAXDIRS 10240
#define MAXFILES 102400

/*
Checks for numeric values in the strings, so that strings like
'image1.img' 'image10.img' 'image11.img' 'image12.img' 'image2.img' etc
are ordered more sensibly.
*/
static int num_strcmp(char *str1, char *str2)
{
	char *p1, *p2;
	int i1, i2;

	for(p1=str1, p2=str2; *p1 && *p2; p1++, p2++)
	{
		if (*p1>='0' && *p1<='9' && *p2>='0' && *p2<='9')
		{
			i1 = i2 = 0;
			while(*p1>='0' && *p1<='9')
				i1 = i1*10 + (*(p1++) -'0');
			while(*p2>='0' && *p2<='9')
				i2 = i2*10 + (*(p2++) -'0');
			if (i1>i2)
				return(1);
			else if (i1<i2)
				return(-1);
		}
		if (*p1 > *p2)
			return(1);
		else if (*p1 < *p2)
			return(-1);
	}
	return(0);
}


static void slowsort(int argc,char **argv)
{
	int i, j;
	for(i=argc; i>1; i--)
		for(j=1; j<i; j++)
		{
			if(num_strcmp(argv[j], argv[j-1]) < 0)
			{
				char *tempp;
				tempp = argv[j];
				argv[j] = argv[j-1];
				argv[j-1] = tempp;
			}
		}
}

static int wildcard(char *attempt, char *actual)
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

static void list2mat(int m, int n, char *list[], mxArray **ptr)
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


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int ndirs = 0, nfiles = 0, len, uid, maxdlen = 0, maxflen = 0;
	long int ngids, gids[128];
	DIR *dirp;
	struct dirent *dp;
	char *filenames[MAXFILES], *directories[MAXDIRS], *ptr, *buf = (char *)0, *bufp, *fullpathname, *filter;
	static struct stat stbuf;

	if ((nrhs != 2) || (nlhs != 2))
		mexErrMsgTxt("Incorrect Usage.");

	if (mxIsNumeric(prhs[0]) || mxIsNumeric(prhs[1]))
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
		/* Apart from automounted directories, the sum of all the lengths
		   of filenames within a directory should not exceed the the size
		   of a directory (as ascertained with stat).  However, for
		   automount directories such as /home, then the size of the
		   directory is equal to the number of files in the directory.
		   This was the reason for a bug that took ages to find.  As a
		   hack, I will add 65536 to the size of buf.  This should protect
		   against most automount directories.  The alternatives would either
		   be a two pass approach to listing the contents of a directory, or
		   to calloc space for files as needed.
		*/
		buf = (char *)mxCalloc((int)stbuf.st_size+65536, 1);
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
						if (ndirs == MAXDIRS)
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
						if (nfiles == MAXFILES)
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

