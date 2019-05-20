/* $Id: shoot_boundary.h 7593 2019-05-20 18:58:16Z john $ */
/* (c) John Ashburner (2011) */

#define BOUND_CIRCULANT 0
#define BOUND_NEUMANN   1
#define BOUND_DIRICHLET 2 
#define BOUND_SLIDING   3

extern mwSignedIndex (*bound)();
extern void set_bound(int t);
extern int  get_bound();

#define BOUND(a,b) bound(a,b)

