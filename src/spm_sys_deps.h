/* $Id: spm_sys_deps.h 112 2005-05-04 18:20:52Z john $ */

#ifndef SYS_DEP
#define SYS_DEP
#ifdef SPM_WIN32
#define rint(x) floor((x)+0.5)  /* round to nearest int */
#define finite(x) mxIsFinite(x) /* finite */
#define SEPCHAR      '\\'          /* directory separator */
#include <process.h>
#else /* SPM_WIN32 */
#define SEPCHAR      '/'           /* directory separator */
extern double rint(double);
#endif /* SPM_WIN32 */
#ifndef MAXNAMLEN
/* default (very long) maximum filename length */
#define MAXNAMLEN 1024
#endif
#endif /* SYS_DEP */
