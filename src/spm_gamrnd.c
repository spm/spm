/*
 * $Id: spm_gamrnd.c 3251 2009-07-06 17:29:44Z guillaume $
 * Guillaume Flandin
 */

/*
 * (Very) Highly inspired from Lightspeed Matlab toolbox, by Tom Minka 
 * http://research.microsoft.com/en-us/um/people/minka/software/lightspeed/
 */

#include <math.h>
#include "mex.h"


/* Return a sample from Uniform on the unit interval
 */
static long ix = 101;
static long iy = 1001;
static long iz = 10001;
double Rand(void)
{
  static float u;
  
  ix = 171*(ix % 177)-2*(ix/177);
  iy = 172*(iy % 176)-2*(iy/176);
  iz = 170*(iz % 178)-2*(iz/178);
  
  if (ix<0) ix = ix + 30269;
  if (iy<0) iy = iy + 30307;
  if (iz<0) iz = iz + 30323;
  
  u = ((float) ix)/30269 +
                ((float) iy)/30307 + ((float) iz)/30323;
  u -= (float)(int)u;
  return(u);
}

/* Returns a sample from Normal(0,1)
 */
double RandN(void)
{
  double x,y,radius;
  /* Generate a random point inside the unit circle */
  do {
    x = 2*Rand()-1;
    y = 2*Rand()-1;
    radius = (x*x)+(y*y);
  } while((radius >= 1.0) || (radius == 0.0));
  /* Box-Muller formula */
  radius = sqrt(-2*log(radius)/radius);
  x *= radius;
  y *= radius;
  return x;
}

/* Returns a sample from Gamma(a, 1).
 */
double GammaRand(double a)
{
  /* Algorithm:
   * G. Marsaglia and W.W. Tsang, A simple method for generating gamma
   * variables, ACM Transactions on Mathematical Software, Vol. 26, No. 3,
   * Pages 363-372, September, 2000.
   * http://portal.acm.org/citation.cfm?id=358414
   */
  double boost, d, c, v;
  if(a < 1) {
    /* boost using Marsaglia's (1961) method: gam(a) = gam(a+1)*U^(1/a) */
    boost = exp(log(Rand())/a);
    a++;
  } 
  else boost = 1;
  d = a-1.0/3; c = 1.0/sqrt(9*d);
  while(1) {
    double x,u;
    do {
      x = RandN();
      v = 1+c*x;
    } while(v <= 0);
    v = v*v*v;
    x = x*x;
    u = Rand();
    if((u < 1-.0331*x*x) || 
       (log(u) < 0.5*x + d*(1-v+log(v)))) break;
  }
  return( boost*d*v );
}

/* Main gateway 
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
  int ndims, len, i;
  int *dims = NULL;
  double a, b;
  double *o = NULL;

  if (nlhs > 1)
    mexErrMsgTxt("Too many output arguments.");
  
  if (nrhs < 2)
    mexErrMsgTxt("Requires at least two input arguments.");

  ndims = (nrhs == 2) ? 1 : nrhs-2;
  dims = (int*)mxMalloc(ndims*sizeof(int));
  
  len = 1;
  for (i=0;i<ndims;i++) {
    dims[i] = (nrhs == 2) ? 1 : (int)mxGetScalar(prhs[i+2]);
    len *= dims[i];
  }
  
  a = mxGetScalar(prhs[0]);
  b = mxGetScalar(prhs[1]);
  
  plhs[0] = mxCreateNumericArray(ndims, dims, mxDOUBLE_CLASS, mxREAL);
  o = mxGetPr(plhs[0]);

  for(i=0;i<len;i++)
    *o++ = b * GammaRand(a);
  
  mxFree(dims);
}
