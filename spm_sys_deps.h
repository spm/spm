/* %W% John Ashburner - from Matthew Brett's work %E% */

#ifndef SYS_DEP
	#define SYS_DEP
	#ifdef SPM_WIN32
		#define rint(x) floor((x)+0.5)  /* round to nearest int */
		#define finite(x) mxIsFinite(x) /* finite */
		#define SEPS      "\\"          /* directory separator */
		#include <process.h>
	#else /* SPM_WIN32 */
		#define SEPS      "/"           /* directory separator */
		extern double rint(double);
	#endif /* SPM_WIN32 */
	#ifndef MAXNAMLEN
		/* default (very long) maximum filename length */
		#define MAXNAMLEN 1024
	#endif
#endif /* SYS_DEP */

