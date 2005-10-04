/* 
 * $Id: spm_make_lookup.h 247 2005-10-04 17:20:34Z guillaume $
 */

/* Generate a lookup table for Lagrange interpolation */

#ifndef _SPM_MAKE_LOOKUP_H_
#define _SPM_MAKE_LOOKUP_H_

void make_lookup_poly(double coord, int q, int dim, int *d1,
	double *table, double **ptpend);
	
void make_lookup_poly_grad(double coord, int q, int dim, int *d1,
	double *table, double *dtable, double **ptpend);
	
void make_lookup_sinc(double coord, int q, int dim, int *d1,
	double *table, double **ptpend);
	
void make_lookup_sinc_grad(double coord, int q, int dim, int *d1,
	double *table, double *dtable, double **ptpend);

#endif /* _SPM_MAKE_LOOKUP_H_ */
