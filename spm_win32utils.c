/* spm_win32utils */
#ifndef lint
static char sccsid[]="@(#)spm_win32utils.c	2.3 Matthew Brett 99/05/17";
#endif

#include <windows.h>
#include <string.h>
#include <dirent.h>
#include "mex.h"

#define STRSIZE 1024

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int len;
	char* actionstr;

	if ((nrhs < 1))
		mexErrMsgTxt("Incorrect Usage.");

	if (mxIsNumeric(prhs[0]))
		mexErrMsgTxt("Action argument must be a string.");

	if (mxGetM(prhs[0]) != 1 )
		mexErrMsgTxt("Action argument must have one row.");

	len = mxGetN(prhs[0]);
	actionstr = (char *)mxCalloc(len+1, sizeof(char));
	mxGetString(prhs[0],actionstr,len+1);

	if (strcmp(actionstr, "username")==0){
		/* get windows 95/8 logon name - does not work for NT */
		DWORD datasize = STRSIZE;
		char regname[STRSIZE];
		DWORD vtdum;
		HKEY regkey;

		RegOpenKeyEx(
 			HKEY_LOCAL_MACHINE,
 			"Network\\Logon",
 		 	0,   // reserved
 		 	KEY_EXECUTE, // security access mask
  			&regkey    // address of handle to open key
		);
		RegQueryValueEx(
 			regkey,
 			"username",
  			NULL,
  			&vtdum,         // address of buffer for value type
  			regname,        // address of data buffer
  			&datasize       // address of data buffer size
		);
		RegCloseKey(
 			 regkey   // handle to key to close
		);

		/* Return registry username */
		plhs[0]=mxCreateString(regname);

	} else if (strcmp(actionstr, "drives")==0){
		int drive, driveno = 0, curdrive;
		UINT prevemode;
		char drives[27];
		char curdir[_MAX_PATH];

		/* turn off insert disk in drive messages */
		prevemode = SetErrorMode(SEM_FAILCRITICALERRORS);

		/* Save current drive. */
		curdrive = _getdrive();
		/* Get the current working directory: */
		_getcwd( curdir, _MAX_PATH );

		/* If we can switch to the drive, it exists. */
		/* start at c:, don't check floppy drives */
		for( drive = 3; drive <= 26; drive++ )
			if( !_chdrive( drive ) )
				drives[driveno++]='A'+drive-1;
		drives[driveno]=0;

		/* Restore original drive, dir and error mode.*/
		_chdrive( curdrive );
	        _chdir( curdir ); 
		SetErrorMode(prevemode);

		/* Return string of drive letters */
		plhs[0]=mxCreateString(drives);
	}
}
