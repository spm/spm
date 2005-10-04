/*
 * $Id: spm_sys_deps.h 247 2005-10-04 17:20:34Z guillaume $
 */

#ifndef _SPM_SYS_DEP_H_
#define _SPM_SYS_DEP_H_

#ifdef SPM_WIN32
  #define SEPCHAR      '\\'  /* directory separator */
#else /* SPM_WIN32 */
  #define SEPCHAR      '/'   /* directory separator */
#endif /* SPM_WIN32 */

#ifndef MAXNAMLEN
#define MAXNAMLEN 1024       /* default (very long) maximum filename length */
#endif

#endif /* _SPM_SYS_DEP_H_ */
