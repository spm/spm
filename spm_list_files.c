#ifndef lint
static char sccsid[]="%W% (c) John Ashburner %E%";
#endif

#include <string.h>
#include <sys/types.h>
#include <sys/param.h>
#include <sys/stat.h>
#include <dirent.h>
#include <unistd.h>

#include "mex.h"

#define MAXDIRS 10240
#define MAXFILES 102400

#include "spm_sys_deps.h"

static int getmask(struct stat *stbuf)
{
	int mask = 0007;
#ifndef SPM_WIN32
        static uid_t uid = -1;
        static gid_t gids[128];
        static int ngids;
	int g;

	if (uid == -1)
	{
		uid   = getuid();
		ngids = getgroups(128,gids);
	}
	if (stbuf->st_uid == uid)
		mask = 0700; /* user */
	else
		for(g=0; g<ngids; g++)
			if (gids[g] == stbuf->st_gid)
				mask = 0070; /* group */
#endif
	return(mask);
}

/* constants for memory alloc to file name storage */
#define MEMBLOCK 65536
#define MAXBUFS 256

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
	for (p1=attempt, p2=actual; ((*p1 == *p2) || (*p1 == '?')) && ((*p1 != '\0') && (*p2 != '\0')); p1++, p2++);
	if ((*p1 == *p2) || ((*p1 == '?') && (*p2 != '\0'))) return 1;
	if  (*p1 != '*') return 0;
	for (;*p1 == '*';p1++);
	for (;(wildcard(p1, p2) == 0) && (*p2 != '\0');p2++);
	if ((*p1 == *p2) || ((*p1 == '?') && (*p2 != '\0'))) return 1;
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
	int ndirs = 0, nfiles = 0, len, maxdlen = 0, maxflen = 0, i;
	DIR *dirp;
	struct dirent *dp;
	char *filenames[MAXFILES], *directories[MAXDIRS], *ptr, *fullpathname, *filter;
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

	if ((stat(fullpathname, &stbuf) != -1) && (dirp = opendir(fullpathname)))
	{
		ptr = fullpathname + strlen(fullpathname);
		(void)strcat(fullpathname, SEPS);
		ptr++;
		for (dp = readdir(dirp); dp != NULL; dp = readdir(dirp))
		{
			*ptr = 0;
			(void)strcpy(ptr, dp->d_name);

			if (stat(fullpathname, &stbuf) != -1
				&& ((stbuf.st_mode & S_IFMT) == S_IFDIR
				|| ((stbuf.st_mode & S_IFMT) == S_IFREG
				&& wildcard(filter,dp->d_name))))
			{
				int mask = getmask(&stbuf);

				if ((stbuf.st_mode & S_IFMT) == S_IFDIR
					&& (mask & 0555 & stbuf.st_mode))
				{
					if (ndirs == MAXDIRS)
						mexErrMsgTxt("Too many subdirectories.");

					len = strlen(dp->d_name);
					if (len > maxdlen) maxdlen = len;
					directories[ndirs] = (char *)mxCalloc(len+1, sizeof(char));
					(void)strcpy(directories[ndirs], dp->d_name);
					ndirs++;
				}
				else if ((stbuf.st_mode & S_IFMT) == S_IFREG
					&&  (mask & 0444 & stbuf.st_mode))
				{
					if (nfiles == MAXFILES)
						mexErrMsgTxt("Too many files match.");

					len = strlen(dp->d_name);
					if (len > maxflen) maxflen = len;
					filenames[nfiles] = (char *)mxCalloc(len+1, sizeof(char));
					(void)strcpy(filenames[nfiles], dp->d_name);
					nfiles++;
				}
			} /* entry statable */
		} /* for entries in directory */
		(void)closedir(dirp);
	} /* directory statable */
		
	slowsort(nfiles, filenames);
	list2mat(nfiles,maxflen,filenames,&plhs[0]);
	
	slowsort(ndirs, directories);
	list2mat(ndirs, maxdlen,directories,&plhs[1]);

	for (i=0;i<nfiles;i++) (void)mxFree((char *)filenames[i]); 
	for (i=0;i<ndirs ;i++) (void)mxFree((char *)directories[i]); 
	(void)mxFree((char *)fullpathname);
}
