#ifndef lint
static char sccsid[] = "%W% (c) John Ashburner MRCCU/FIL %E%";
#endif lint

#include "cmex.h"
#include "volume.h"
#include <math.h>
extern double floor(), fabs();


/*
INPUTS
T[3*nz*ny*nx + ni*4] - current transform

dt2           - data type of image to normalize
scale2        - scalefactor of image to normalize
dim2[3]       - x, y & z dimensions of image to normalize
dat2[dim2[2]*dim2[1]*dim2[0]]     - voxels of image to normalize

ni            - number of templates
dt1[ni]       - data types of templates
scale1[ni]    - scalefactors of templates
dim1[3]       - x, y & z dimensions of templates
dat1[ni*dim1[2]*dim1[1]*dim1[0]] - voxels of templates

nx            - number of basis functions in x
BX[dim1[0]*nx]- basis functions in x
ny            - number of basis functions in y 
BY[dim1[1]*ny]- basis functions in y
nz            - number of basis functions in z
BZ[dim1[2]*nz]- basis functions in z

M[4*4]        - transformation matrix
samp[3]       - frequency of sampling template.


OUTPUTS
alpha[(3*nz*ny*nx+ni*4)^2] - A'A
 beta[(3*nz*ny*nx+ni*4)]   - A'b
pss[1*1]                   - sum of squares difference
pnsamp[1*1]                - number of voxels sampled

*/


void mrqcof(T, alpha, beta, pss, dt2,scale2,dim2,dat2, ni,dt1,scale1,dim1,dat1, nx,BX, ny,BY, nz,BZ, M, samp, pnsamp)
double T[], alpha[], beta[], pss[], BX[], BY[], BZ[], M[], scale2, scale1[];
int dt2, dim2[], dt1[], dim1[], ni, nx, ny, nz, samp[], *pnsamp;
unsigned char *dat1[], dat2[];
{
	int i1,i2, s0[3], x1,x2, y1,y2, z1,z2, m1, m2, nsamp = 0, ni4;
	double dvds0[3],dvds1[3], *dvdt, s2[3], *ptr1, *ptr2, *Tz, *Ty, tmp,dtmp,
		*betaxy, *betax, *alphaxy, *alphax, ss=0.0, *scale1a;

	ni4 = ni*4;

	/* rate of change of voxel with respect to change in parameters */
	dvdt    = (double *)mxCalloc( 3*nx    + ni4                 , sizeof(double));

	/* Intermediate storage arrays */
	Tz      = (double *)mxCalloc( 2*3*nx*ny                     , sizeof(double));
	Ty      = (double *)mxCalloc( 2*3*nx                        , sizeof(double));
	betax   = (double *)mxCalloc( 3*nx    + ni4                 , sizeof(double));
	betaxy  = (double *)mxCalloc( 3*nx*ny + ni4                 , sizeof(double));
	alphax  = (double *)mxCalloc((3*nx    + ni4)*(3*nx    + ni4), sizeof(double));
	alphaxy = (double *)mxCalloc((3*nx*ny + ni4)*(3*nx*ny + ni4), sizeof(double));

	/* pointer to scales for each of the template images */
	scale1a = T + 3*nx*ny*nz;

	/* only zero half the matrix */
	m1 = 3*nx*ny*nz+ni4;
	for (x1=0;x1<m1;x1++)
	{
		for (x2=0;x2<=x1;x2++)
			alpha[m1*x1+x2] = 0.0;
		beta[x1]= 0.0;
	}

	for(s0[2]=0; s0[2]<dim1[2]-1; s0[2]+=samp[2]) /* For each plane of the template images */
	{
		/* build up the deformation field from it's seperable form */
		for(i1=0, ptr1=T, ptr2=Tz; i1<3; i1++, ptr1 += nz*ny*nx, ptr2 += ny*nx)
			for(x1=0; x1<nx*ny; x1++)
			{
				dtmp = tmp = 0.0;
				for(z1=0; z1<nz; z1++)
				{
					tmp  += ptr1[x1+z1*ny*nx] *  BZ[dim1[2]*z1+s0[2]];
					if (i1==2)
						dtmp += ptr1[x1+z1*ny*nx] * BZ[dim1[2]*z1+s0[2]+1];					
				}

				ptr2[x1] = tmp;

				if (i1==2)
					ptr2[x1 + 3*nx*ny] = dtmp;
				else
					ptr2[x1 + 3*nx*ny] = tmp;
			}

		/* only zero half the matrix */
		m1 = 3*nx*ny+ni4;
		for (x1=0;x1<m1;x1++)
		{
			for (x2=0;x2<=x1;x2++)
				alphaxy[m1*x1+x2] = 0.0;
			betaxy[x1]= 0.0;
		}

		for(s0[1]=0; s0[1]<dim1[1]-1; s0[1]+=samp[1]) /* For each row of the template images plane */
		{
			/* build up the deformation field (and gradients) from it's seperable form */
			for(i1=0, ptr1=Tz, ptr2=Ty; i1<3; i1++, ptr1+=ny*nx, ptr2+=nx)
			{
				for(x1=0; x1<nx; x1++)
				{
					dtmp = tmp = 0.0;
					for(y1=0; y1<ny; y1++)
					{
						tmp  += ptr1[x1+y1*nx] *  BY[dim1[1]*y1+s0[1]];
						if (i1==1)
							dtmp += ptr1[x1+nx*(y1+3*ny)] * BY[dim1[1]*y1+s0[1]+1];
						else
							dtmp += ptr1[x1+nx*(y1+3*ny)] * BY[dim1[1]*y1+s0[1]  ];

					}
					ptr2[x1     ] = tmp;
					ptr2[x1+3*nx] = dtmp;
				}
			}

			/* only zero half the matrix */
			m1 = 3*nx+ni4;
			for(x1=0;x1<m1;x1++)
			{
				for (x2=0;x2<=x1;x2++)
					alphax[m1*x1+x2] = 0.0;
				betax[x1]= 0.0;
			}

			for(s0[0]=0; s0[0]<dim1[0]-1; s0[0]+=samp[0]) /* For each pixel in the row */
			{
				double trans[3], dtrans[3];

				/* nonlinear deformation of the template space, followed by the affine transform */
				for(i1=0, ptr1 = Ty; i1<3; i1++, ptr1 += nx)
				{
					dtmp = tmp = 1.0;
					for(x1=0; x1<nx; x1++)
					{
						tmp  += ptr1[x1] * BX[dim1[0]*x1+s0[0]];
						if (i1==0)
							dtmp += ptr1[x1+3*nx] * BX[dim1[0]*x1+s0[0]+1];
						else
							dtmp += ptr1[x1+3*nx] * BX[dim1[0]*x1+s0[0]  ];
					}
					 trans[i1] = tmp + s0[i1];
					dtrans[i1] = 1 + dtmp - tmp; /* derivative */
				}

				/* Affine component */
				s2[0] = M[0+4*0]*trans[0] + M[0+4*1]*trans[1] + M[0+4*2]*trans[2] + M[0+4*3];
				s2[1] = M[1+4*0]*trans[0] + M[1+4*1]*trans[1] + M[1+4*2]*trans[2] + M[1+4*3];
				s2[2] = M[2+4*0]*trans[0] + M[2+4*1]*trans[1] + M[2+4*2]*trans[2] + M[2+4*3];


				/* is the transformed position in range? */
				if (s2[0]>=1.0 && s2[0]<dim2[0] && s2[1]>=1.0 &&
					s2[1]<dim2[1] && s2[2]>=1.0 && s2[2]<dim2[2])
				{
					double val000, val001, val010, val011, val100, val101, val110, val111;
					double dx1, dy1, dz1, dx2, dy2, dz2;
					double v,dv;
					int xcoord, ycoord, zcoord, off1, off2;

					nsamp ++;

					/* coordinates to resample from */
					xcoord = (int)floor(s2[0]);
					ycoord = (int)floor(s2[1]);
					zcoord = (int)floor(s2[2]);

					/* Calculate the interpolation weighting factors*/
					dx1=s2[0]-xcoord; dx2=1.0-dx1;
					dy1=s2[1]-ycoord; dy2=1.0-dy1;
					dz1=s2[2]-zcoord; dz2=1.0-dz1;

					/* get pixel values of 8 nearest neighbours */
					off1 = xcoord-1 + dim2[0]*(ycoord-1  + dim2[1]*(zcoord-1));
					switch (dt2)
					{
					case UNSIGNED_CHAR:
						val111 = dat2[off1]; val011 = dat2[off1+1];
						off2 = off1+dim2[0];
						val101 = dat2[off2]; val001 = dat2[off2+1];
						off1+= dim2[0]*dim2[1];
						val110 = dat2[off1]; val010 = dat2[off1+1];
						off2 = off1+dim2[0];
						val100 = dat2[off2]; val000 = dat2[off2+1];
						break;
					case SIGNED_SHORT:
						val111 = ((short *)dat2)[off1]; val011 = ((short *)dat2)[off1+1];
						off2 = off1+dim2[0];
						val101 = ((short *)dat2)[off2]; val001 = ((short *)dat2)[off2+1];
						off1+= dim2[0]*dim2[1];
						val110 = ((short *)dat2)[off1]; val010 = ((short *)dat2)[off1+1];
						off2 = off1+dim2[0];
						val100 = ((short *)dat2)[off2]; val000 = ((short *)dat2)[off2+1];
						break;
					case SIGNED_INT:
						val111 = ((int *)dat2)[off1]; val011 = ((int *)dat2)[off1+1];
						off2 = off1+dim2[0];
						val101 = ((int *)dat2)[off2]; val001 = ((int *)dat2)[off2+1];
						off1+= dim2[0]*dim2[1];
						val110 = ((int *)dat2)[off1]; val010 = ((int *)dat2)[off1+1];
						off2 = off1+dim2[0];
						val100 = ((int *)dat2)[off2]; val000 = ((int *)dat2)[off2+1];
						break;
					case FLOAT:
						val111 = ((float *)dat2)[off1]; val011 = ((float *)dat2)[off1+1];
						off2 = off1+dim2[0];
						val101 = ((float *)dat2)[off2]; val001 = ((float *)dat2)[off2+1];
						off1+= dim2[0]*dim2[1];
						val110 = ((float *)dat2)[off1]; val010 = ((float *)dat2)[off1+1];
						off2 = off1+dim2[0];
						val100 = ((float *)dat2)[off2]; val000 = ((float *)dat2)[off2+1];
						break;
					case DOUBLE:
						val111 = ((double *)dat2)[off1]; val011 = ((double *)dat2)[off1+1];
						off2 = off1+dim2[0];
						val101 = ((double *)dat2)[off2]; val001 = ((double *)dat2)[off2+1];
						off1+= dim2[0]*dim2[1];
						val110 = ((double *)dat2)[off1]; val010 = ((double *)dat2)[off1+1];
						off2 = off1+dim2[0];
						val100 = ((double *)dat2)[off2]; val000 = ((double *)dat2)[off2+1];
						break;
					default:
						mexErrMsgTxt("Bad data type.");
					}

					/* resampled pixel value */
					v =     (((val111*dx2 + val011*dx1)*dy2
						+ (val101*dx2 + val001*dx1)*dy1)*dz2
						+((val110*dx2 + val010*dx1)*dy2
						+ (val100*dx2 + val000*dx1)*dy1)*dz1)*scale2;

					/* local gradients accross resampled pixel (in space of object image) */
					dvds0[0] = ((dy2*(val011-val111) + dy1*(val001-val101))*dz2
						 + (dy2*(val010-val110) + dy1*(val000-val100))*dz1)*scale2;

					dvds0[1] = ((dx2*(val101-val111) + dx1*(val001-val011))*dz2
						 + (dx2*(val100-val110) + dx1*(val000-val010))*dz1)*scale2;

					dvds0[2] = ((dx2*(val110-val111) + dx1*(val010-val011))*dy2
						 + (dx2*(val100-val101) + dx1*(val000-val001))*dy1)*scale2;

					/* transform the gradients to the same space as the template */
					dvds1[0] = (M[0+4*0]*dvds0[0] + M[1+4*0]*dvds0[1] + M[2+4*0]*dvds0[2])*dtrans[0];
					dvds1[1] = (M[0+4*1]*dvds0[0] + M[1+4*1]*dvds0[1] + M[2+4*1]*dvds0[2])*dtrans[1];
					dvds1[2] = (M[0+4*2]*dvds0[0] + M[1+4*2]*dvds0[1] + M[2+4*2]*dvds0[2])*dtrans[2];


					/* there is no change in the contribution from BY and BZ, so only
					   work from BX */
					for(i1=0; i1<3; i1++)
					{
						for(x1=0; x1<nx; x1++)
							dvdt[i1*nx+x1] = -dvds1[i1] * BX[dim1[0]*x1+s0[0]];
					}

					/* coordinate of template voxel to sample */
					off2 = s0[0] + dim1[0]*(s0[1] + dim1[1]*s0[2]);

					/* dv will be the difference between the object image
					   and the linear combination of templates */
					dv = v;
					for(i1=0; i1<ni; i1++)
					{
						double tmp;
						switch (dt1[i1])
						{
						case UNSIGNED_CHAR:
							tmp = ((unsigned char *)(dat1[i1]))[off2]*scale1[i1];
							break;
						case SIGNED_SHORT:
							tmp = ((short         *)(dat1[i1]))[off2]*scale1[i1];
							break;
						case SIGNED_INT:
							tmp = ((int           *)(dat1[i1]))[off2]*scale1[i1];
							break;
						case FLOAT:
							tmp = ((float         *)(dat1[i1]))[off2]*scale1[i1];
							break;
						case DOUBLE:
							tmp = ((double        *)(dat1[i1]))[off2]*scale1[i1];
							break;
						default:
							mexErrMsgTxt("Bad data type.");
						}

						/* linear combination of image and image modulated by constant
						   gradients in x, y and z */
						dvdt[i1*4  +3*nx] = tmp;
						dvdt[i1*4+1+3*nx] = tmp*s2[0];
						dvdt[i1*4+2+3*nx] = tmp*s2[1];
						dvdt[i1*4+3+3*nx] = tmp*s2[2];

						dv -= dvdt[i1*4  +3*nx]*scale1a[i1*4];
						dv -= dvdt[i1*4+1+3*nx]*scale1a[i1*4+1];
						dv -= dvdt[i1*4+2+3*nx]*scale1a[i1*4+2];
						dv -= dvdt[i1*4+3+3*nx]*scale1a[i1*4+3];
					}


					/* cf Numerical Recipies "mrqcof.c" routine */
					m1 = 3*nx+ni4;
					for(x1=0; x1<m1; x1++)
					{
						for (x2=0;x2<=x1;x2++)
							alphax[m1*x1+x2] += dvdt[x1]*dvdt[x2];
						betax[x1] += dvdt[x1]*dv;
					}

					/* sum of squares */
					ss += dv*dv;
				}
			}

			m1 = 3*nx*ny+ni4;
			m2 = 3*nx+ni4;

			/* Kronecker tensor products */
			for(y1=0; y1<ny; y1++)
			{
				double wt = BY[dim1[1]*y1+s0[1]];

				for(i1=0; i1<3; i1++)	/* loop over deformations in x, y and z */
				{
					/* spatial-spatial covariances */
					for(i2=0; i2<=i1; i2++)	/* symmetric matrixes - so only work on half */
					{
						for(y2=0; y2<=y1; y2++)
						{
							/* Kronecker tensor products with BY'*BY */
							double wt2 = wt * BY[dim1[1]*y2+s0[1]];

							ptr1 = alphaxy + nx*(m1*(ny*i1 + y1) + ny*i2 + y2);
							ptr2 = alphax  + nx*(m2*i1 + i2);

							for(x1=0; x1<nx; x1++)
							{
								for (x2=0;x2<=x1;x2++)
									ptr1[m1*x1+x2] += wt2 * ptr2[m2*x1+x2];
							}
						}
					}

					/* spatial-intensity covariances */
					ptr1 = alphaxy + nx*(m1*ny*3 + ny*i1 + y1);
					ptr2 = alphax  + nx*(m2*3 + i1);
					for(x1=0; x1<ni4; x1++)
					{
						for (x2=0;x2<nx;x2++)
							ptr1[m1*x1+x2] += wt * ptr2[m2*x1+x2];
					}

					/* spatial component of beta */
					for(x1=0; x1<nx; x1++)
						betaxy[x1+nx*(ny*i1 + y1)] += wt * betax[x1 + nx*i1];
				}
			}
			ptr1 = alphaxy + nx*(m1*ny*3 + ny*3);
			ptr2 = alphax  + nx*(m2*3 + 3);
			for(x1=0; x1<ni4; x1++)
			{
				/* intensity-intensity covariances  */
				for (x2=0; x2<=x1; x2++)
					ptr1[m1*x1 + x2] += ptr2[m2*x1 + x2];

				/* intensity component of beta */
				betaxy[nx*ny*3 + x1] += betax[nx*3 + x1];
			}
		}

		m1 = 3*nx*ny*nz+ni4;
		m2 = 3*nx*ny+ni4;

		/* Kronecker tensor products */
		for(z1=0; z1<nz; z1++)
		{
			double wt = BZ[dim1[2]*z1+s0[2]];

			for(i1=0; i1<3; i1++)	/* loop over deformations in x, y and z */
			{
				/* spatial-spatial covariances */
				for(i2=0; i2<=i1; i2++)	/* symmetric matrixes - so only work on half */
				{
					for(z2=0; z2<=z1; z2++)
					{
						/* Kronecker tensor products with BZ'*BZ */
						double wt2 = wt * BZ[dim1[2]*z2+s0[2]];

						ptr1 = alpha   + nx*ny*(m1*(nz*i1 + z1) + nz*i2 + z2);
						ptr2 = alphaxy + nx*ny*(m2*i1 + i2);
						for(y1=0; y1<ny*nx; y1++)
						{
							for (y2=0;y2<=y1;y2++)
								ptr1[m1*y1+y2] += wt2 * ptr2[m2*y1+y2];
						}
					}
				}

				/* spatial-intensity covariances */
				ptr1 = alpha   + nx*ny*(m1*nz*3 + nz*i1 + z1);
				ptr2 = alphaxy + nx*ny*(m2*3 + i1);
				for(y1=0; y1<ni4; y1++)
				{
					for (y2=0;y2<ny*nx;y2++)
						ptr1[m1*y1+y2] += wt * ptr2[m2*y1+y2];
				}

				/* spatial component of beta */
				for(y1=0; y1<ny*nx; y1++)
					beta[y1 + nx*ny*(nz*i1 + z1)] += wt * betaxy[y1 + nx*ny*i1];
			}
		}

		ptr1 = alpha   + nx*ny*(m1*nz*3 + nz*3);
		ptr2 = alphaxy + nx*ny*(m2*3 + 3);
		for(y1=0; y1<ni4; y1++)
		{
			/* intensity-intensity covariances */
			for(y2=0;y2<=y1;y2++)
				ptr1[m1*y1 + y2] += ptr2[m2*y1 + y2];

			/* intensity component of beta */
			beta[nx*ny*nz*3 + y1] += betaxy[nx*ny*3 + y1];
		}

		mexPrintf(".");

	}


	/* Fill in the symmetric bits
	   - OK I know some bits are done more than once - but it shouldn't matter. */

	m1 = 3*nx*ny*nz+ni4;
	for(i1=0; i1<3; i1++)	
	{
		double *ptrz, *ptry, *ptrx;
		for(i2=0; i2<=i1; i2++)
		{
			ptrz = alpha + nx*ny*nz*(m1*i1 + i2);
			for(z1=0; z1<nz; z1++)
				for(z2=0; z2<=z1; z2++)
				{
					ptry = ptrz + nx*ny*(m1*z1 + z2);
					for(y1=0; y1<ny; y1++)
						for (y2=0;y2<=y1;y2++)
						{
							ptrx = ptry + nx*(m1*y1 + y2);
							for(x1=0; x1<nx; x1++)
								for(x2=0; x2<x1; x2++)
									ptrx[m1*x2+x1] = ptrx[m1*x1+x2];
						}
					for(x1=0; x1<nx*ny; x1++)
						for (x2=0; x2<x1; x2++)
							ptry[m1*x2+x1] = ptry[m1*x1+x2];
				}
			for(x1=0; x1<nx*ny*nz; x1++)
				for (x2=0; x2<x1; x2++)
					ptrz[m1*x2+x1] = ptrz[m1*x1+x2];
		}
	}
	for(x1=0; x1<nx*ny*nz*3; x1++)
		for (x2=0; x2<x1; x2++)
			alpha[m1*x2+x1] = alpha[m1*x1+x2];

	ptr1 = alpha + nx*ny*nz*3;
	ptr2 = alpha + nx*ny*nz*3*m1;
	for(x1=0; x1<nx*ny*nz*3; x1++)
		for(x2=0; x2<ni4; x2++)
			ptr1[m1*x1+x2] = ptr2[m1*x2+x1];

	ptr1 = alpha + nx*ny*nz*3*(m1+1);
	for(x1=0; x1<ni4; x1++)
		for (x2=0; x2<x1; x2++)
			ptr1[m1*x2+x1] = ptr1[m1*x1+x2];

	*pss = ss;
	*pnsamp = nsamp;

	mxFree((char *)dvdt);
	mxFree((char *)Tz);
	mxFree((char *)Ty);
	mxFree((char *)betax);
	mxFree((char *)betaxy);
	mxFree((char *)alphax);
	mxFree((char *)alphaxy);

	mexPrintf("\n");
}



#ifdef __STDC__
void mexFunction(int nlhs, Matrix *plhs[], int nrhs, Matrix *prhs[])
#else
mexFunction(nlhs, plhs, nrhs, prhs)
int nlhs, nrhs;
Matrix *plhs[], *prhs[];
#endif
{
	extern double rint();
	MAPTYPE **map1, *map2;
	int i, nx,ny,nz,ni=1, samp[3];
	int dim1[3], dim2[3], dt1[32], dt2, nsamp;
        double *M, *BX, *BY, *BZ, *T, scale1[32], scale2, fwhm;
	unsigned char *dat1[32], *dat2;

        if (nrhs != 8 || nlhs > 3)
        {
                mexErrMsgTxt("Inappropriate usage. ([A,B,var]=f(V1,V2,M,BX,BY,BZ,T,fwhm);)");
        }

        map1 = get_maps(prhs[0], &ni);
        map2 = get_map(prhs[1]);

	dim1[0]   = map1[0]->xdim;
	dim1[1]   = map1[0]->ydim;
	dim1[2]   = map1[0]->zdim;


	for (i=0; i<8; i++)
		if (!mxIsNumeric(prhs[i]) || mxIsComplex(prhs[i]) ||
			!mxIsFull(prhs[i]) || !mxIsDouble(prhs[i]))
			mexErrMsgTxt("Inputs must be numeric, real, full and double.");

	if (mxGetM(prhs[2]) != 4 || mxGetN(prhs[2]) != 4)
		mexErrMsgTxt("Transformation matrix must be 4x4.");
	M = mxGetPr(prhs[2]);


	if ( mxGetM (prhs[3]) != dim1[0]) mexErrMsgTxt("Wrong sized X basis functions.");
	nx = mxGetN (prhs[3]);
	BX = mxGetPr(prhs[3]);

	if ( mxGetM (prhs[4]) != dim1[1]) mexErrMsgTxt("Wrong sized Y basis functions.");
	ny = mxGetN (prhs[4]);
	BY = mxGetPr(prhs[4]);

	if ( mxGetM (prhs[5]) != dim1[2]) mexErrMsgTxt("Wrong sized Z basis functions.");
	nz = mxGetN (prhs[5]);
	BZ = mxGetPr(prhs[5]);

	T = mxGetPr(prhs[6]);
	if (mxGetM(prhs[6])*mxGetN(prhs[6]) != 3*nx*ny*nz+ni*4)
		mexErrMsgTxt("Transform is wrong size.");

	if ( mxGetM(prhs[7])*mxGetN(prhs[7]) != 1) mexErrMsgTxt("FWHM should be a single value.");
	fwhm = mxGetPr(prhs[7])[0];

	/* sample about every fwhm/2 */
	samp[0]   = rint(fwhm/2.0/map1[0]->xpixdim); samp[0] = ((samp[0]<1) ? 1 : samp[0]);
	samp[1]   = rint(fwhm/2.0/map1[0]->ypixdim); samp[1] = ((samp[1]<1) ? 1 : samp[1]);
	samp[2]   = rint(fwhm/2.0/map1[0]->zpixdim); samp[2] = ((samp[2]<1) ? 1 : samp[2]);

	for(i=0; i<ni; i++)
	{
		dat1[i]   = map1[i]->data;
		dt1[i]    = map1[i]->datatype;
		scale1[i] = map1[i]->scalefactor;
		if (dim1[0] != map1[i]->xdim || dim1[1] != map1[i]->ydim || dim1[2] != map1[i]->zdim)
			mexErrMsgTxt("Volumes must have same dimensions.");
	}

	dat2      = map2->data;
	dim2[0]   = map2->xdim;
	dim2[1]   = map2->ydim;
	dim2[2]   = map2->zdim;
	dt2       = map2->datatype;
	scale2    = map2->scalefactor;


        plhs[0] = mxCreateFull(3*nx*ny*nz+ni*4,3*nx*ny*nz+ni*4,REAL);
        plhs[1] = mxCreateFull(3*nx*ny*nz+ni*4,1,REAL);
	plhs[2] = mxCreateFull(1,1,REAL);

	mrqcof(T, mxGetPr(plhs[0]), mxGetPr(plhs[1]), mxGetPr(plhs[2]),
		dt2,scale2,dim2,dat2,
		ni,dt1,scale1,dim1,dat1,
		nx, BX, ny, BY, nz, BZ, M, samp, &nsamp);

	mxGetPr(plhs[2])[0] = mxGetPr(plhs[2])[0] / (nsamp - (3*nx*ny*nz + ni*4));
	free_maps(map1);
}

