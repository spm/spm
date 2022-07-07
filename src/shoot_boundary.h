/*
 * John Ashburner
 * Copyright (C) 2011-2022 Wellcome Centre for Human Neuroimaging
 */

#define BOUND_CIRCULANT 0
#define BOUND_NEUMANN   1
#define BOUND_DIRICHLET 2 
#define BOUND_SLIDING   3

extern mwSignedIndex (*bound)();
extern void set_bound(int t);
extern int  get_bound();

#define BOUND(a,b) bound(a,b)
