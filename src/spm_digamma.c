/*
 * $Id: spm_digamma.c 112 2005-05-04 18:20:52Z john $
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mex.h"

/*
 *  Mex file derived from a FORTRAN program by D. E. Amos,
 *  ACM Transactions on Mathematical Software, 1983.
 *  Obtained from NETLIB: http://www.netlib.org/toms/610.
 */

#ifndef MAX
#define MAX(A,B) (((A)>(B))?(A):(B))
#endif

#ifndef MIN
#define MIN(A,B) (((A)<(B))?(A):(B))
#endif

void dpsifn(double x, int n, int kode, int m, double *ans) {
	int i, j, k, mx, nmax, nn, np, nx;
	double arg, den, eps, fln, fn, fnp, fns, fx, rln, rxsq, r1m4, r1m5, s, slope, t, ta, tk, tol, tols, trm[22], trmr[100], tss, tst, tt, t1, t2, wdtol, xdmln, xdmy, xinc, xln, xm, xmin, xq, yint, xinf;
	/* Bernoulli numbers */
    double b[] = { 1.00000000000000000e+00, -5.00000000000000000e-01, \
                   1.66666666666666667e-01, -3.33333333333333333e-02, \
                   2.38095238095238095e-02, -3.33333333333333333e-02, \
                   7.57575757575757576e-02, -2.53113553113553114e-01, \
                   1.16666666666666667e+00, -7.09215686274509804e+00, \
                   5.49711779448621554e+01, -5.29124242424242424e+02, \
                   6.19212318840579710e+03, -8.65802531135531136e+04, \
                   1.42551716666666667e+06, -2.72982310678160920e+07, \
                   6.01580873900642368e+08, -1.51163157670921569e+10, \
                   4.29614643061166667e+11, -1.37116552050883328e+13, \
                   4.88332318973593167e+14, -1.92965793419400681e+16 };
				   
	nmax = 100;
	nx = 0; /* to remove a warning stating that `nx' might be used uninitialized... */
	xinf = mxGetInf();
	
	if (x < 0)
		mexErrMsgTxt("X must be nonnegative.");
	if (n < 0)
		mexErrMsgTxt("Internal error (N)");
	if (kode != 1)
		mexErrMsgTxt("Internal error (KODE)");
	if (m < 1)
		mexErrMsgTxt("Internal error (M)");
		
	if (x == 0.0) {
		for (i=0;i<m;i++) ans[i] = xinf;
		return;
	}
	else if (mxIsInf(x)) {
		for (i=0;i<m;i++) ans[i] = xinf;
		if (n == 0) ans[0] = -xinf;
		return;
	}
	else if (mxIsNaN(x)) {
		for (i=0;i<m;i++) ans[i] = x;
		return;
	}
	
	nn = n + m - 1;
	fn = (double)nn;
	fnp = fn + 1.0;
	r1m5 = log10(2.0);
	r1m4 = 1.110223e-16;
	wdtol = MAX(r1m4, 0.5e-18);
	xln = log(x);
	if (x < wdtol) goto line260;
	
	/* Compute xmin and the number of terms of the series, fln+1 */
	rln = r1m5 * 53.0;
	rln = MIN(rln, 18.06);
	fln = MAX(rln, 3.0) - 3.0;
	yint = 3.5 + 0.4 * fln;
	slope = 0.21 + fln * (0.0006038 * fln + 0.008677);
	xm = yint + slope * fn;
	mx = (int)xm + 1;
	xmin = (double)mx;
	if (n != 0) {
		xm = -2.302 * rln - MIN(0.0,xln);
		fns = (double)n;
		arg = xm / fns;
		arg = MIN(0.0,arg);
		eps = exp(arg);
		xm = 1.0 - eps;
		if (fabs(arg) < 1.0e-3) xm = -arg;
		fln = x * xm / eps;
		xm = xmin - x;
		if ((xm > 7.0) && (fln < 15.0)) goto line200;
	}
	xdmy = x;
	xdmln = xln;
	xinc = 0.0;
	if (x <= xmin) {
		nx = (int)x;
		xinc = xmin - (double)nx;
		xdmy = x + xinc;
		xdmln = log(xdmy);
	}
	
	/* Generate w(n+m-1,x) by the asymptotic expansion */
	t = fn * xdmln;
	t1 = xdmln + xdmln;
	t2 = t + xdmln;
	tss = exp(-t);
	tt = 0.5 / xdmy;
	t1 = tt;
	tst = wdtol * tt;
	if (n != 0) t1 = tt + 1.0 / fn;
	rxsq = 1.0 / (xdmy * xdmy);
	ta = 0.5 * rxsq;
	t = fnp * ta;
	s = t * b[2];
	
	if (fabs(s) >= tst) {
		tk = 2.0;
		for (k=3;k<22;k++) {
			t = t * ((tk + fn + 1.0) / (tk + 1.0)) * ((tk + fn) / (tk + 2.0)) * rxsq;
			trm[k] = t * b[k];
			if (fabs(trm[k]) < tst) break;
			s += trm[k];
			tk += 2.0;
		}
	}
	s = (s + t1) * tss;
	if (xinc == 0.0) goto line100;
	
	/* backward recur from xdmy to x */
	nx = (int)xinc;
	np = nn + 1;
	if (nx > nmax) 
		mexErrMsgTxt("Internal error (NMAX)");
	if (nn == 0) goto line160;
	xm = xinc - 1.0;
	fx = x + xm;
	
	/* this loop should not be changed. fx is accurate when x is small */
	for (i=0;i<nx;i++) {
		trmr[i] = pow(fx,-np);
		s += trmr[i];
		xm -= 1.0;
		fx = x + xm;
	}
	line100:
	ans[m-1] = s;
	if (fn == 0.0) goto line180;
	
	/* Generate lower derivatives, j<n+m-1 */
	if (m == 1) return;
	for (j=1;j<m;j++) {
		fnp = fn;
		fn -= 1.0;
		tss *= xdmy;
		t1 = tt;
		if (fn != 0) t1 = tt + 1.0 / fn;
		t = fnp * ta;
		s = t * b[2];
		if (fabs(s) >= tst) {
			tk = 3.0 + fnp;
			for (k=3;k<22;k++) {
				trm[k] = trm[k] * fnp / tk;
				if (fabs(trm[k]) < tst) break;
				s += trm[k];
				tk += 2.0;
			}
		}
		s = (s + t1) * tss;
		if (xinc == 0.0) goto line140;
		if (fn == 0) goto line160;
		xm = xinc - 1.0;
		fx = x + xm;
		for (i=0;i<nx;i++) {
			trmr[i] *= fx;
			s += trmr[i];
			xm -= 1.0;
			fx = x + xm;
		}
		line140:
		mx = m - j - 1;
		ans[mx] = s;
		if (fn == 0.0) goto line180;
	}
	return;
	/* recursion for n = 0 */
	line160:
	for (i=1;i<=nx;i++)
		s += 1.0 / (x + (double)(nx - i));
	line180:
	if (kode != 2) {
		ans[0] = s - xdmln;
		return;
	}
	if (xdmy == x) return;
	xq = xdmy / x;
	ans[0] = s - log(xq);
	return;
	
	/* compute by series (x+k)^(-(n+1)), k=0,1,2,... */
	line200:
	nn = (int)fln + 1;
	np = n + 1;
	t1 = (fns + 1.0) * xln;
	t = exp(-t1);
	s = t;
	den = x;
	for (i=0;i<nn;i++) {
		den += 1.0;
		trm[i] = pow(den,-np);
		s += trm[i];
	}
	ans[0] = s;
	if (n != 0)
		if (kode == 2)
			ans[0] = s + xln;
	if (m == 1) return;
	
	/* Generate higher derivatives, J > n */
	tol = wdtol / 5.0;
	for (j=1;j<m;j++) {
		t /= x;
		s = t;
		tols = t * tol;
		den = x;
		for (i=0;i<nn;i++) {
			den += 1.0;
			trm[i] = trm[i] / den;
			s += trm[i];
			if (trm[i] < tols) break;
		}
		ans[j] = s;
	}
	return;
		
	/* Small x < unit round off */
	line260:
	ans[0] = pow(x,-n-1);
	if (m != 1)
		for (k=0,i=1;i<m;i++,k++)
			ans[k+1] = ans[k] / x;
	if (n == 0) return;
	if (kode == 2) ans[0] += xln;
	return;
}

void unscale(int k0, int m, int n, double *y) {
	int i, j, k;
	double s = -1.0;
	for (k=1;k<=k0;k++)
		s = -k * s;
	k = k0;
	for (i=0;i<m;i++) {
		for (j=0;j<n;j++)
			y[i+m*j] = s * y[i+m*j];
		k++;
		s *= -k;
	}
}

/* --- GATEWAY FUNCTION --- */
void mexFunction(int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[]) {
	
	int i, m, n, k0, k1, mn;
	double *pk = NULL, *px = NULL, *py = NULL;

	if (nrhs < 1)
		mexErrMsgTxt("Not enough input arguments.");
	if (nrhs > 2)
		mexErrMsgTxt("Too many input arguments.");
	if (nrhs == 1) {
		if (mxIsComplex(prhs[0]))
			mexErrMsgTxt("Argument must be real.");
		m = mxGetM(prhs[0]);
		n = mxGetN(prhs[0]);
		plhs[0] = mxCreateDoubleMatrix(m,n,mxREAL);
		px = mxGetPr(prhs[0]);
		py = mxGetPr(plhs[0]);
		mn = m * n;
		for (i=0;i<mn;i++)
			dpsifn(*(px + i), 0, 1, 1, py + i);
		unscale(0, 1, mn, py);
	}
	else {
		if ((mxIsComplex(prhs[0])) || (mxIsComplex(prhs[1])))
			mexErrMsgTxt("Arguments must be real.");
		m = mxGetM(prhs[0]) * mxGetN(prhs[0]);
		n = mxGetM(prhs[1]) * mxGetN(prhs[1]);
		plhs[0] = mxCreateDoubleMatrix(m,n,mxREAL);
		if ((m == 0) || (n == 0)) return;
		pk = mxGetPr(prhs[0]);
		px = mxGetPr(prhs[1]);
		py = mxGetPr(plhs[0]);
		k0 = pk[0];
		k1 = pk[m-1];
		if (k0 < 0)
			mexErrMsgTxt("Order must be nonnegative.");
		if (k1 != k0+m-1)
			mexErrMsgTxt("Order must be consecutive integers.");
		for (i=0;i<n;i++)
			dpsifn(*(px + i), k0, 1, m, py + i * m);
		unscale(k0, m, n, py);
	}
}
