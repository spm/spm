/* spm_win32utils */
#ifndef lint
static char sccsid[]="%W% Matthew Brett %E%";
#endif

#include <windows.h>
#include <string.h>

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
		plhs[0]=mxCreateString(regname);
	}
}
