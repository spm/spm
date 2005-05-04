/* $Id: spm_mapping.h 112 2005-05-04 18:20:52Z john $
   matlab dependent high level data access and map manipulation routines */

#include <sys/types.h>
#ifdef SPM_WIN32
#include <windows.h>
#include <memory.h>
#else
#include <sys/mman.h>
#endif

#include "spm_vol_access.h"
#include "mex.h"

void free_maps(MAPTYPE *maps, int n);

MAPTYPE *get_maps(const mxArray *ptr, int *n);

void voxdim(MAPTYPE *map, double vdim[3]);

int get_dtype(const mxArray *ptr);

