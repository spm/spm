/* $Id: shoot_boundary.h 4583 2011-12-06 16:03:01Z john $ */
/* (c) John Ashburner (2011) */

#define BOUND_CIRCULANT 0
#define BOUND_NEUMANN 1

extern mwSignedIndex (*bound)();
extern void set_bound(int t);
extern int  get_bound();

#define BOUND(a,b) bound(a,b)

