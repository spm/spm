/* spm_mapping.h */

/* matlab dependent high level data access and map
manipulation routines */

#include "spm_vol_access.h"
#include "mex.h"

void free_maps(MAPTYPE *maps, int n);

MAPTYPE *get_maps(const mxArray *ptr, int *n);

void voxdim(MAPTYPE *map, double vdim[3]);

int get_dtype(const mxArray *ptr);

