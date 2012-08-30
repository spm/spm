/* $Id: shoot_boundary.h 4875 2012-08-30 20:04:30Z john $ */
/* (c) John Ashburner (2011) */

#define BOUND_CIRCULANT 0
#define BOUND_NEUMANN 1

extern mwSignedIndex (*bound)();
extern void set_bound(int t);
extern int  get_bound();

#define BOUND(a,b) bound(a,b)

