/*
 * John Ashburner, Mikael Brudfors & Yael Balbastre
 * Copyright (C) 2020-2022 Wellcome Centre for Human Neuroimaging
 */
#include "spm_mex.h"
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#define EXP(x) fastexp(x)

typedef struct
{
    double *s2;
    double *s1;
    double *s0;
} SStype;

typedef struct
{
    mwSize  P;
    double *mu;
    double *b;
    double *W;
    double *nu;
    double *gam;
    double *conN;
    double *conT;
} GMMtype;

typedef struct
{
    mwSize  Po;
    mwSize  Pm;
    unsigned char *observed;
    double *Wt;
    double *L_mm;
    double *m_m;
    double *m_o;
    unsigned char *obs;
} MissInfType;

static const double pi = 3.1415926535897931;
#define MaxChan   ((mwSize)50) /* largest integer valued float is 2^52 */
#define Undefined ((mwSize)0xFFFFFFFFFFFFF)


/* A (hopefully) faster approximation to exp.
 * Note that some precision is lost for values
 * of x further from integers. */
static double fastexp(double x)
{
    double r, rr;
    mwSignedIndex i;
    static double lkp_mem[256], *exp_lkp = lkp_mem+128;

    /* exp(i+r) = exp(i)*exp(r), where:
     *     exp(i) (where i is an integer) is from the lookup table;
     *     exp(r) (residual) is from a generalised continued fraction
     *            https://en.wikipedia.org/wiki/Exponential_function#Continued_fractions_for_ex
     *
     * Should not encounter values more extreme than -128 or 127,
     * particularly as the upper limit of x will be 0 and values
     * of x below log(eps)=-36.04 should be numerically equivalent.*/
    i  = (mwSignedIndex)rint(x);
    if (i<-128) i = -128;
    if (i> 127) i =  127;
    if (exp_lkp[i]==0.0) exp_lkp[i] = exp((double)i);

    r  = x - (double)i;
    rr = r*r;
/*  return exp_lkp[i] * (1.0+2.0*r/(2.0-r+rr/(6.0+rr/(10.0+rr/14.0))));
 *  return exp_lkp[i] * (1.0+2.0*r/(2.0-r+rr/(6.0+rr/(10.0)))); */
    return exp_lkp[i] * (1.0+2.0*r/(2.0-r+rr/6.0));
}

static mwSize is_observed(mwSize code, mwSize i)
{
    return (code>>i) & (mwSize)1;
}

static mwSize num_observed(mwSize code, mwSize P)
{
    mwSize i, Po;
    for(i=0, Po=0; i<P; i++)
        Po += is_observed(code,i);
    return Po;
}


/* Allocate memory for a data structure, and assign pointers from
 * the structure to the appropriate parts of s0_ptr, s1_ptr and s2_ptr.
 * P - number of image channels
 * K - number of Gaussians
 * s0_ptr, s1_ptr & s2_ptr - memory allocated for storing sufficient
 *                           statistics.
 * Returns the allocated data structure, with pointers assigned.
 */
static /*@null@*/ SStype *suffstat_pointers(mwSize P, mwSize K, double *s0_ptr, double *s1_ptr, double *s2_ptr)
{
    SStype /*@NULL@*/ *suffstat;
    suffstat = (SStype *)calloc((size_t)1<<(size_t)P,sizeof(SStype));
    if (suffstat != NULL)
    {
        mwSize code,o0=0,o1=0,o2=0;
        for(code=0; code<((mwSize)1<<P); code++)
        {
            mwSize Po;
            Po = num_observed(code,P);
            suffstat[code].s0 = &(s0_ptr[o0]);
            suffstat[code].s1 = &(s1_ptr[o1]);
            suffstat[code].s2 = &(s2_ptr[o2]);
            o0 += K;
            o1 += K*Po;
            o2 += K*Po*Po;
        }
    }
    return suffstat;
}


/* Read vector of mean and variance.
 * N1 - number of voxels each volume of mf and mv
 * P  - number of volumes in mf and mv
 * mf - E[f],         size N1 x P
 * vf - diag(Var[f]), size N1 x P
 * x  - extracted means,     length P
 * v  - extracted variances, length P
 * Returns a code indicating which volumes have a finite value (ie data not missing)
 */
static mwSize get_vox(mwSize N1, mwSize P, float mf[], float vf[], /*@out@*/ double x[], /*@out@*/ double v[])
{
    mwSize j, j1, o, code;
    for(j=0, j1=0, o=0, code=0; j<P; j++, o+=N1)
    {
        double tmp = (double)mf[o];
        if (isfinite(tmp)!=0)
        {
            x[j1] = tmp;
            v[j1] = vf[o];
            code |= (mwSize)1<<j;
            j1++;
        }
    }
    return code;
}


/* log-sum-exp function: log(sum(exp(q)))
 * K - length of q
 *
static double lse(mwSize K, double q[])
{
    mwSize k;
    double mx, s;
    for(k=1, mx=q[0]; k<K; k++)
    {
        if (q[k]>mx) mx = q[k];
    }
    for(k=0, s=0.0; k<K; k++)
        s += EXP(q[k]-mx);
    return log(s)+mx;
}
 */


/* softmax function: exp(q)/sum(exp(q))
 * K - length of q or p
 * q - input data
 * p - output data
 * Input and output could be the same
 * Returns lse(q)
 */
static double softmax1(mwSize K, double q[], /*@out@*/ double p[])
{
    mwSize k;
    double mx, s;
    for(k=1, mx=q[0]; k<K; k++)
        if (q[k]>mx) mx = q[k];
    for(k=0, s=0.0; k<K; k++)
        s += (p[k] = EXP(q[k]-mx));
    for(k=0; k<K; k++)
        p[k] /= s;
    return log(s) + mx;
}


/* softmax function: exp(q)/(sum(exp(q)) + 1)
 * K - length of q or p
 * q - input data
 * p - output data
 * Input and output could be the same
 * Returns lse([q 0])
 */
static double /*@unused@*/ softmax(mwSize K, double q[], /*@out@*/ double p[])
{
    mwSize k;
    double mx, s;
    for(k=0, mx=0; k<K; k++)
        if (q[k]>mx) mx = q[k];
    for(k=0, s=EXP(-mx); k<K; k++)
        s += (p[k] = EXP(q[k]-mx));
    for(k=0; k<K; k++)
        p[k] /= s;
    return log(s) + mx;
}


/* \Del^2 = (\mu-x)^T W (\mu-x) + trace(W diag(v))
 * P  - dimensions of mu, W, x & v
 * mu - mean, size P x 1
 * W  - covariance, size P x P
 * x  - expectation of data, size P x 1
 * v  - variance of data, size P x 1 (ie diagonal of covariance)
 */
static double del2(mwSize P, double mu[], double W[], double x[], double v[])
{
    mwSize j,i;
    double d=0.0, r, *wj;
    for(j=0,wj=W; j<P; j++, wj+=P)
    {
        r  = x[j]-mu[j];
        d += wj[j]*(r*r+v[j]);
        for(i=j+1; i<P; i++)
            d += 2.0*r*wj[i]*(x[i]-mu[i]);
    }
    return d;
}


/* psi / digamma function
 * From http://web.science.mq.edu.au/~mjohnson/code/digamma.c
 * by Mark Johnson, Professor of Language Sciences (CORE),
 * Department of Computing, Faculty of Science,
 * Macquarie University, Sydney, Australia. 
 */
static double psi(double z)
{
    double f = 0, r, r2, r4;
    /* psi(z) = psi(z+1) - 1/z */
    for (f=0; z<7.0; z++) f -= 1.0/z;

    z -= 1.0/2.0;
    r  = 1.0/z;
    r2 = r*r;
    r4 = r2*r2;
    f += log(z)+(1.0/24.0)*r2-(7.0/960.0)*r4+(31.0/8064.0)*r4*r2-(127.0/30720.0)*r4*r4;
    return f;
}

/* Add log-likelihoods from Normal distributions, accounting for uncertainty of
 * the parameters
 *
 * K    - number of Gaussians.
 * gmm  - data structure for Gaussian distributions with different
 *        combinations of missing data.
 * code - Indicates which data are missing.
 * x    - E[x]
 * v    - Var[x]
 * p    - on input: a vector of logs of prior probabilities
 *        on output: logs of likelihoods are added to the log priors
 */
static void Nloglikelihoods(mwSize K, GMMtype gmm[], mwSize code, double x[], double v[], double p[])
{
    mwSize P, k;
    double *mu, *b, *W, *nu, *gam, *con;
    P   = gmm[code].P;
    mu  = gmm[code].mu;
    b   = gmm[code].b;
    W   = gmm[code].W;
    nu  = gmm[code].nu;
    gam = gmm[code].gam;
    con = gmm[code].conN;

    for(k=0; k<K; k++, W+=P*P, mu+=P)
        p[k] += con[k] - 0.5*nu[k]*del2(P, mu, W, x, v);
}

/* Add log-likelihoods from Student's T distributions, accounting
 * for uncertainty of the parameters
 *
 * K    - number of Gaussians.
 * gmm  - data structure for Gaussian distributions with different
 *        combinations of missing data.
 * code - Indicates which data are missing.
 * x    - E[x]
 * v    - Var[x]
 * p    - on input: a vector of logs of prior probabilities
 *        on output: logs of likelihoods are added to the log priors
 */
static void Tloglikelihoods(mwSize K, GMMtype gmm[], mwSize code, double x[], double v[], double p[])
{
    mwSize P, k;
    double *mu, *b, *W, *nu, *con;
    /*
       Compute other log-likelihoods from a mixture of Student's t distributions.
       See Eqns. 10.78-10.82 & B.68-B.72 in Bishop's PRML book.
       In practice, it only improves probabilities by a tiny amount.

       ln St(x|mu,Lam,tau)
       gammaln((tau + P)/2.0) - gammaln(tau/2.0) + ld1/2 - (P/2.0)*log(tau*pi) - ((tau+P)/2)*log(1 + del2/tau)
       where:
           Lam  = W*(nu+1-P)*beta/(1+beta)
           tau  = nu+1-P
           del2 = (x-mu)'*Lam*(x-mu)
           ld1  = log(det(Lam)) = log(det(W)) + P*log((nu+1-P)*beta/(1+beta))
       This gives:
       gammaln((nu+1)/2.0) - gammaln((nu+1-P)/2.0) + (ld + P*log((nu+1-P)*beta/(1+beta)))/2 - ...
       (P/2.0)*log((nu+1-P)*pi) - ...
       ((nu+1)/2)*log(1 + beta/(beta+1)*((x-mu)'*W*(x-mu) + sum(diag(W).*vf)))
    */
    P   = gmm[code].P;
    mu  = gmm[code].mu;
    b   = gmm[code].b;
    W   = gmm[code].W;
    nu  = gmm[code].nu;
    con = gmm[code].conT;

    for(k=0; k<K; k++, W+=P*P, mu+=P)
        p[k] += con[k] - 0.5*(nu[k]+1.0)*log(1.0 + b[k]/(b[k]+1.0)*del2(P, mu, W, x, v));
}


/*
 * lnP = E[\ln A] 
 * E[\ln {\bf a}_k] = \psi({\bf a}_k) - \psi(\sum_i a_{ik})
 *
 */
static void Dloglikelihoods(mwSize i, mwSize K, double *lnP, double *p)
{
    mwSize k;
    for(k=0; k<K; k++)
        p[k] += lnP[i+256*k];
}

/* Construct a vector of log tissue priors for a voxel
 *
 * N1  - Number of voxels in each 3D volume.
 * lp  - Pointer to first voxel in the volumes.
 * K   - Number of Gaussian distributions
 * lkp - Lookup table indicating which Gaussian is
 *       associated with which tissue class.
 * p   - Output vector of log tissue priors.
 */
static int logpriors(mwSize N1, float *lp, mwSize K, mwSize *lkp, /*@out@*/ double *p)
{
    mwSize k;
/*  double l; */
    for(k=0; k<K; k++)
    {
        double lpk;
        lpk = (double)lp[N1*lkp[k]];
        if (isfinite(lpk)==0)
            return 0;
        p[k] = lpk;
    }
    /*
    l = lse(K,p);
    for(k=0; k<K; k++)
        p[k] -= l;
    */
    return 1;
}


/* Cholesky decomposition
 * n  - dimension of matrix a
 * a  - an n \times n matrix
 * p  - an n \times 1 vector
 *
 * A triangle of the input matrix is partially overwritten
 * by the output. Diagonal elements are stored in p.
 */
static void choldc(mwSize n, double a[], /*@out@*/ double p[])
{
    mwSize    i, j;
    mwSignedIndex k;
    double sm, sm0;

    sm0  = 1e-40;
    for(i=0; i<n; i++) sm0 = sm0 + a[i*n+i];
    sm0 *= 1e-7;
    sm0 *= sm0;

    for(i=0; i<n; i++)
    {
        for(j=i; j<n; j++)
        {
            sm = a[i*n+j];
            for(k=(mwSignedIndex)i-1; k>=0; k--)
               sm -= a[i*n+k] * a[j*n+k];
            if(i==j)
            {
                if(sm <= sm0) sm = sm0;
                p[i] = sqrt(sm);
            }
            else
                a[j*n+i] = sm / p[i];
        }
    }
}


/* Solve a least squares problem with the results from a
 * Cholesky decomposition
 *
 * n     - Dimension of matrix and data.
 * a & p - Cholesky decomposed matrix.
 * b     - Vector of input data.
 * x     - Vector or outputs.
 */
static void cholls(mwSize n, const double a[], const double p[],
            const double b[], /*@out@*/ double x[])
{
    mwSignedIndex i, k;
    double sm;

    for(i=0; i<(mwSignedIndex)n; i++)
    {
        sm = b[i];
        for(k=i-1; k>=0; k--)
            sm -= a[i*n+k]*x[k];
        x[i] = sm/p[i];
    }
    for(i=(mwSignedIndex)n-1; i>=0; i--)
    {
        sm = x[i];
        for(k=i+1; k<(mwSignedIndex)n; k++)
            sm -= a[k*n+i]*x[k];
        x[i] = sm/p[i];
    }
}


/* n! */
static mwSize factorial(mwSize n)
{
    static mwSize products[21];
    if (products[0]==0)
    {
        mwSize i;
        products[0] = 1;
        for(i=1; i<21; i++)
            products[i] = products[i-1]*i;
    }
    return products[n];
}


/* Compute space required for storing sufficient statistics.
 *
 * P              - Number of image volumes.
 * K              - Number of tissue classes.
 * *m0, *m1 & *m2 - Space needed for the zeroeth,
 *                  first and second moments.
 */
void space_needed(mwSize P, mwSize K, mwSize *m0, mwSize *m1, mwSize *m2)
{
    mwSize m;
    for(m=0, *m0=0, *m1=0, *m2=0; m<=P; m++)
    {
        mwSize nel;
        nel = K*factorial(P)/(factorial(m)*factorial(P - m));
        *m0 += nel;
        *m1 += nel*m;
        *m2 += nel*m*m;
    }
}

/* Allocate memory for a data structure for representing
 * GMMs with missing data
 *
 * P - Number of images/channels.
 * K - Number of Gaussians
 */
static /*@null@*/ GMMtype *allocate_gmm(mwSize P, mwSize K)
{
    mwSize o, code, i, n0=0,n1=0,n2=0;
    double *buf;
    unsigned char *bytes;
    GMMtype /*@NULL@*/ *gmm;
    space_needed(P, K, &n0, &n1, &n2);

    o     = ((mwSize)1<<P)*sizeof(GMMtype);
    bytes = calloc((size_t)o+(size_t)(n0*(mwSize)5+n1+n2)*sizeof(double),1);
    gmm   = (GMMtype *)bytes;
    if (gmm!=NULL)
    {
        buf   = (double *)(bytes + o);
        o     = 0;
        for(code=0; code<((mwSize)1<<P); code++)
        {
            mwSize nel = 0;
            for(i=0; i<code; i++) nel += is_observed(code,i);
            gmm[code].P    = nel;
            gmm[code].mu   = buf+o; o += K*nel;
            gmm[code].b    = buf+o; o += K;
            gmm[code].W    = buf+o; o += K*nel*nel;
            gmm[code].nu   = buf+o; o += K;
            gmm[code].gam  = buf+o; o += K;
            gmm[code].conN = buf+o; o += K;
            gmm[code].conT = buf+o; o += K;
        }
    }
    return gmm;
}

/* Invert a matrix
 *
 * P - Matrix dimensions
 * W - Matrix (input, P \times P)
 * S - Matrix inverse (output, P \times P)
 * T - Scratch space (P*(P+1))
 */
static double invert(mwSize P, double *W /* P*P */, double *S /* P*P */, double *T /* P*(P+1) */)
{
    mwSize i, j, PP=P*P;
    double ld = 0.0, *p;
    for(i=0; i<PP; i++) T[i] = W[i];
    p = T+PP;
    choldc(P,T,p);
    for(j=0; j<P; j++)
    {
       ld += log(p[j]);
        /* Column of identity matrix */
        for(i=0; i<P; i++) S[i+j*P]=0.0;
        S[j+j*P] = 1.0;

        cholls(P, T, p, S+j*P, S+j*P);
    }
    return -2.0*ld;
}

static /*@null@*/ MissInfType *allocate_missinf(mwSize P, mwSize K)
{
    mwSize o, code, m, mem;
    unsigned char *bytes;
    MissInfType /*@NULL@*/ *missinf;

    for(m=0, mem=0; m<=P; m++)
    {
        mwSize nel;
        nel = factorial(P)/(factorial(m)*factorial(P - m));
        mem += nel*K*(P-m)*m*sizeof(double); /* L_mo */
        mem += nel*K*P*sizeof(double);
        mem += nel*K*m*m*sizeof(double);
        mem += nel*((P+sizeof(double)-1)/sizeof(double))*sizeof(double);
    }

    o     = ((mwSize)1<<P)*sizeof(MissInfType);
    bytes = calloc((size_t)(o + mem), (size_t)1);

    missinf = (MissInfType *)bytes;
    if (missinf!=NULL)
    {
        double *buf;
 
        buf   = (double *)(bytes + o);
        o     = 0;
        for(code=0; code<((mwSize)1<<P); code++)
        {
            mwSize Po = num_observed(code,P);
            mwSize Pm = P-Po;
            missinf[code].Po   = Po;
            missinf[code].Pm   = Pm;
            missinf[code].Wt   = buf+o; o += K*Po*Pm;
            missinf[code].L_mm = buf+o; o += K*Pm*Pm;
            missinf[code].m_m  = buf+o; o += K*Pm;
            missinf[code].m_o  = buf+o; o += K*Po;
            missinf[code].obs  = (unsigned char *)(buf+o);
            o += (P+sizeof(double)-1)/sizeof(double);
        }
    }
    return missinf;
}


static void dispmat(double *f, mwSize d1, mwSize d2, mwSize d3)
{
    mwSize i, j, k;
    for(k=0; k<d3; k++)
    {
        for(i=0; i<d1; i++)
        {
            for(j=0; j<d2; j++)
                printf(" %g", f[i+d1*(d2*k+j)]);
            printf("\n");
        }
        printf("--\n");
    }
}


static MissInfType /*@NULL@*/ *prepare_missinf(mwSize P, mwSize K, double *L /* P*P */, double *m)
{
    MissInfType *missinf;
    mwSize code;

    missinf = allocate_missinf(P, K);
    if (missinf==NULL) return missinf;

    for(code=0; code<(mwSize)1<<P; code++)
    {
        mwSize i, j, im, io, jm, jo, k;
        mwSize Po, Pm;
        double p[64], *Lmm, *Lmo, *mm, *mo;
        unsigned char *obs;

        Po  = missinf[code].Po;
        Pm  = missinf[code].Pm;
        Lmm = missinf[code].L_mm; /* partially decomposed on output */
        mm  = missinf[code].m_m;
        mo  = missinf[code].m_o;
        Lmo = missinf[code].Wt;
        obs = missinf[code].obs;

        for(j=0, jm=0, jo=0; j<P; j++)
        {
            if (is_observed(code,j)==0)
            {/* missing */
                obs[j] = (unsigned char)0;
                jm ++;
            }
            else
            {/* observed */
                obs[j] = (unsigned char)1;
                jo ++;
            }
        }

        for(k=0; k<K; k++)
        {
            for(i=0, im=0, io=0; i<P; i++)
            {
                if (is_observed(code,i)==0)
                {
                    for(j=0, jo=0,jm=0; j<P; j++)
                    {
                        if (is_observed(code,j)==0)
                            Lmm[im+Pm*(jm++)] = L[i+P*(j+P*k)];
                        else
                            Lmo[im+Pm*(jo++)] = L[i+P*(j+P*k)];
                    }
                    mm[im++] = m[i+P*k];
                }
                else
                    mo[io++] = m[i+P*k];
            }

            /* Wt = Lmm\Lmo */
            choldc(Pm,Lmm,p);
            for(jo=0; jo<Po; jo++)
                cholls(Pm,Lmm, p, Lmo+Pm*jo, Lmo+Pm*jo);

            /* Update pointers to k'th vector/matrix */
            mm  += Pm;
            mo  += Po;
            Lmm += Pm*Pm;
            Lmo += Pm*Po;
        }
    }

/* For debugging 
    for(code=0; code<(mwSize)1<<P; code++)
    {
        mwSize Po, Pm;
        Po  = missinf[code].Po;
        Pm  = missinf[code].Pm;
        printf("\ncode=%d Po=%d Pm=%d\n", code, Po, Pm);
        printf("\nLmm\n"); dispmat(missinf[code].L_mm,Pm,Pm,K);
        printf("\nmm\n");  dispmat(missinf[code].m_m,1,Pm,K);
        printf("\nmo\n");  dispmat(missinf[code].m_o,1,Po,K);
        printf("\nLmo\n"); dispmat(missinf[code].Wt,Pm,Po,K);
    }
*/
    return missinf;
}


/* Construct a data structure for storing a variational Gaussian
 * mixture model for handling missing data.
 *
 * P  - Dimension
 * K  - Number of Gaussians
 * mu,b,W,nu - Variational Bayesian GMM parameters
 *             mu - P \times K
 *             b  - 1 \times K
 *             W  - P \times P \times K
 *             nu - 1 \times K
 * gam       - Mixing proportions (1 \times K).
 *
 * The function returns the data structure.
 */
static /*@null@*/ GMMtype *sub_gmm(mwSize P, mwSize K, double *mu, double *b, double *W, double *nu, double *gam)
{
    const double log2pi = log(2*pi), log2 = log(2.0);
    double *S, *Si;
    GMMtype *gmm;
    mwSize k, code, PP = P*P;

    if ((gmm = allocate_gmm(P,K)) == NULL) return gmm;
    if ((S   = (double *)calloc((size_t)(P*((mwSize)3*P+(mwSize)1)),sizeof(double))) == NULL)
    {
       (void)free((void *)gmm);
       return NULL;
    }
    Si = S + PP;
    for(k=0; k<K; k++)
    {
        double lgam = log(gam[k]);
        (void)invert(P,W+PP*k,S,S+PP);
        for(code=0; code<(mwSize)1<<P; code++)
        {
            mwSize j,j1, Po;
            double ld, ld1, eld;
            Po               = gmm[code].P;
            gmm[code].nu[k]  = nu[k] - (double)(P-Po);
            gmm[code].b[k]   = b[k];
            gmm[code].gam[k] = lgam;
            for(j=0, j1=0; j<P; j++)
            {
                if ((((mwSize)1<<j) & code) != 0)
                {
                    mwSize i, i1;
                    gmm[code].mu[j1+Po*k] = mu[j+P*k];
                    for(i=0, i1=0; i<P; i++)
                    {
                        if ((((mwSize)1<<i) & code) != 0)
                        {
                            Si[i1+Po*j1] = S[i+P*j];
                            i1++;
                        }
                    }
                    j1++;
                }
            }
            ld = invert(Po,Si,gmm[code].W+k*Po*Po,Si+Po*Po);

            /* Constant term for VB mixture of Gaussians
               E[ln N(x | m, L^{-1})] w.r.t. Gaussian-Wishart */
            for(j=0,eld=0.0; j<Po; j++) eld += psi((gmm[code].nu[k]-(double)j)*0.5);
            eld    += Po*log2 + ld;
            gmm[code].conN[k] = 0.5*(eld - (double)Po*(log2pi+1.0/gmm[code].b[k])) + lgam;

            /* Constant term for VB mixture of T distributions */
            ld1 = ld + (double)Po*log((gmm[code].nu[k]+1.0-(double)Po)*b[k]/(b[k]+1.0));
            gmm[code].conT[k] = lgamma(0.5*(gmm[code].nu[k]+1.0)) - lgamma(0.5*(gmm[code].nu[k]+1.0-(double)Po)) +
                                0.5*ld1 - 0.5*Po*log((gmm[code].nu[k]+1.0-(double)Po)*pi) + lgam;
        }
    }
    (void)free((void *)S);
    return gmm;
}

/* Compute sufficient statistics in a way that handles missing data
 *
 * nf   - Vector of dimensions (n_x, n_y, n_z, P).
 * mf   - E[f], dimensions nf.
 * vf   - Var[f], dimensions nf.
 * gmm  - Gaussian mixture model data structure.
 * nm   - Dimensions of log tissue priors (4 elements).
 * skip - Sampling density of tissue priors (in x, y and z).
 * lkp  - Lookup table relating Gaussians to tissue classes.
 * lp   - Log tissue priors
 * suffstat - Data stucture to hold resulting sufficient statistics.
 */
static double suffstats_missing(mwSize nf[], float mf[], float vf[], unsigned char label[],
                      mwSize K, GMMtype gmm[], double lnP[],
                      mwSize nm[], mwSize skip[], mwSize lkp[], float lp[],
                      SStype suffstat[], double H[])
{
    mwSize K1, i0,i1,i2, n2,n1,n0, P, Nf, Nm, code;
    double ll = 0.0, mx[MaxChan], vx[MaxChan], p[128];

    P  = nf[3];
    Nf = nf[0]*nf[1]*nf[2];
    K1 = nm[3];
    Nm = nm[0]*nm[1]*nm[2];

    n2 = nm[2]/skip[2]; if (n2>nf[2]) n2 = nf[2];
    n1 = nm[1]/skip[1]; if (n1>nf[1]) n1 = nf[1];
    n0 = nm[0]/skip[0]; if (n0>nf[0]) n0 = nf[0];

    for(i2=0; i2<n2; i2++)
    {
        for(i1=0; i1<n1; i1++)
        {
            mwSize off_f, off_m;
            off_f = nf[0]*(i1         + nf[1]*i2);
            off_m = nm[0]*(i1*skip[1] + nm[1]*i2*skip[2]);
            for(i0=0; i0<n0; i0++)
            {
                mwSize i, im;
                i    = i0         + off_f;
                im   = i0*skip[0] + off_m;
                code = get_vox(Nf,P,mf+i,vf+i,mx,vx);
                if (code>0 && logpriors(Nm, lp+im, K, lkp, p)!=0)
                {
                    mwSize j, j1, k, Po;
                    double *s0, *s1, *s2;
                    Nloglikelihoods(K, gmm, code, mx, vx, p);

                    if (label!=NULL)
                    {
                        mwSize labi = (mwSize)(label[im]);
                        Dloglikelihoods(labi, K, lnP, p);
                        ll += softmax1(K,p,p);
                        for(k=0; k<K; k++)
                            H[labi + 256*k] += p[k];
                    }
                    else
                        ll += softmax1(K,p,p);

                    Po  = gmm[code].P;
                    s0  = suffstat[code].s0;
                    s1  = suffstat[code].s1;
                    s2  = suffstat[code].s2;
                    for(k=0; k<K; k++, s2+=Po*Po, s1+=Po, s0++)
                    {
                        double pk = p[k];
                        *s0 += pk;
                        for(j=0; j<Po; j++)
                        {
                            double mxj = mx[j];
                            double px  = pk*mxj;
                            s1[j]      += px;
                            s2[j+Po*j] += pk*(mxj*mxj+vx[j]);
                            for(j1=j+1; j1<Po; j1++)
                                s2[j1+Po*j] += px*mx[j1];
                        }
                    }
                }
            }
        }
    }

    /* Add in upper triangle second order sufficiant statistics */
    for(code=1; code<((mwSize)1<<P); code++)
    {
        mwSize j, j1, k, Po;
        double *s2;
        Po = gmm[code].P;
        s2 = suffstat[code].s2;
        for(k=0; k<K; k++, s2+=Po*Po)
        {
            for(j=0; j<Po; j++)
            {
                for(j1=j+1; j1<Po; j1++)
                    s2[j+Po*j1] = s2[j1+Po*j];
            }
        }
    }
    return ll;
}


/* Constructs a gmm structure from mu, b, W, nu, and gam, as well as 
 * a structure containing pointers to the sufficient statistics.
 * It then calls suffstats_missing before freeing up the structures.
 */
double call_suffstats_missing(mwSize nf[], float mf[], float vf[], unsigned char label[],
    mwSize K, double mu[], double b[], double W[], double nu[], double gam[], double lnP[],
    mwSize nm[], mwSize skip[], mwSize lkp[], float lp[],
    double s0_ptr[], double s1_ptr[], double s2_ptr[], double H[])
{
    mwSize P = nf[3];
    GMMtype *gmm;
    SStype  *suffstat;
    double ll=0.0;

    if (P>=MaxChan || K>=128) return NAN;

    if ((gmm      = sub_gmm(P, K, mu, b, W, nu, gam))==NULL) return NAN;
    if ((suffstat = suffstat_pointers(P, K, s0_ptr, s1_ptr, s2_ptr)) == NULL)
    {
        (void)free((void *)gmm);
        return NAN;
    }
    ll = suffstats_missing(nf, mf, vf, label, K, gmm, lnP, nm, skip, lkp, lp, suffstat, H);
    (void)free((void *)gmm);
    (void)free((void *)suffstat);
    return ll;
}


/* Compute responsibilities in a way that handles missing data.
 * Responsibilities used for fitting the GMM are constructed
 * from a VB GMM, whereas those not used ar constructed from
 * a VB mixture of Student's T distributions.
 *
 * nf   - Vector of dimensions (n_x, n_y, n_z, P).
 * mf   - E[f], dimensions nf.
 * vf   - Var[f], dimensions nf.
 * gmm  - Gaussian mixture model data structure.
 * nm   - Dimensions of log tissue priors (4 elements).
 * skip - Sampling density for GMM vs TMM (in x, y and z).
 * lkp  - Lookup table relating Gaussians to tissue classes.
 * lp   - Log tissue priors
 * r    - Responsibilities (n_x, n_y, n_z, max(lkp)).
 */
static double responsibilities(mwSize nf[], mwSize skip[], float mf[], float vf[], unsigned char label[],
              mwSize K, GMMtype *gmm, double lnP[],
              mwSize K1, mwSize lkp[], float lp[],
              float r[])
{
    mwSize P, N1, i0,i1,i2;
    double ll = 0.0, mx[MaxChan], vx[MaxChan], p[128];

    P  = nf[3];
    N1 = nf[0]*nf[1]*nf[2];

    for(i2=0; i2<nf[2]; i2++)
    {
        for(i1=0; i1<nf[1]; i1++)
        {
            mwSize off_f;
            off_f = nf[0]*(i1 + nf[1]*i2);
            for(i0=0; i0<nf[0]; i0++)
            {
                mwSize i, code, k, k1;
                i    = i0+off_f;
                code = get_vox(N1,P,mf+i,vf+i,mx,vx);
                if (logpriors(N1, lp+i, K, lkp, p)!=0)
                {
                    if (code!=0)
                    {
                        if (label!=NULL) Dloglikelihoods((mwSize)(label[i]), K, lnP, p);
                        if ((i2%skip[2])==0 && ((i1%skip[1])==0) & ((i0%skip[0])==0))
                            Nloglikelihoods(K, gmm, code, mx, vx, p);
                        else
                            Tloglikelihoods(K, gmm, code, mx, vx, p);
                        ll += softmax1(K,p,p);
                        for(k=0; k<K; k++)
                        {
                            k1 = lkp[k];
                            if (k1<K1-1)
                                r[i+k1*N1] += p[k];
                        }
                    }
                    else
                    {
                        for(k1=0; k1<K1-1; k1++)
                            r[i+k1*N1]  = NAN;
                        /*
                        (void)softmax(K1,p,p);
                        for(k1=0; k1<K1-1; k1++)
                            r[i+k1*N1] += p[k1];
                         */
                    }
                }
                else
                    for(k1=0; k1<K1-1; k1++)
                        r[i+k1*N1]  = NAN;
            }
        }
    }
    return ll;
}


/* Constructs a gmm structure from mu, b, W, nu, and gam, using this
 * to call responsibilities.
 */
double call_responsibilities(mwSize nf[], mwSize skip[], float mf[], float vf[], unsigned char label[],
    mwSize K, double mu[], double b[], double W[], double nu[], double gam[], double lnP[],
    mwSize K1, mwSize lkp[], float lp[],
    float r[])
{
    mwSize P = nf[3];
    GMMtype *gmm;
    double ll;

    if (P>=MaxChan || K>=128) return NAN;
    if ((gmm      = sub_gmm(P, K, mu, b, W, nu, gam))==NULL) return NAN;

    ll = responsibilities(nf, skip, mf, vf, label, K, gmm, lnP, K1, lkp, lp, r);

    (void)free((void *)gmm);
    return ll;
}

/* Gradient and Hessian for INU updates
 *
 * The computations (two channels only) can be checked with

% Some MATLAB Symbolic Toolbox working...
syms w_11 w_12 w_22 mu_1 mu_2 x_1 x_2 b_1 b_2 mx_1 mx_2 real
syms vx_1 vx_2 positive
W  = [w_11 w_12; w_12 w_22]; % Precision of Gaussian
mu = [mu_1; mu_2];           % Mean of Gaussian
x  = [x_1; x_2];
mx = [mx_1; mx_2];           % E[x]
B  = diag([b_1; 0]);         % INU as a funciton of b_1

% Objective function for a single Gaussian. Extending to more is trivial.
E0  = (x-expm(-B)*mu)'*(expm(B)'*W*expm(B))*(x-expm(-B)*mu)/2 - log(det(expm(B)'*W*expm(B)))/2;

% The above objective function is equivalent to:
E = (expm(B)*x-mu)'*W*(expm(B)*x-mu)/2 - log(det(expm(B)'*W*expm(B)))/2;

% We're using a VB approach with x ~ N(mx,diag(vx)), so compute the expected E.
pdf1 = sym('1/sqrt(2*pi*vx_1)*exp(-(x_1-mx_1)^2/(2*vx_1))');          % x_1 ~ N(mx_1,vx_1)
pdf2 = sym('1/sqrt(2*pi*vx_2)*exp(-(x_2-mx_2)^2/(2*vx_2))');          % x_2 ~ N(mx_2,vx_2)
E    = simplify(int(int(E*pdf1*pdf2,x_1,-Inf,Inf),x_2,-Inf,Inf),1000) % Expectation (takes a while)

% We now assume mx is the expectation of the INU corrected image according to the
% old parameters and vx is the expected variance.  Because
% exp(b+b_old)*x = exp(b)*exp(b_old)*x, we can now assume our initial estimates for
% b are zero and treat exp(b_old)*x as x.
% A quadratic approximation (around b=0) is obtained by:
E0     = simplify(subs(E,b_1,0),1000);
G0     = simplify(subs(diff(E,b_1),b_1,0),1000);           % Gradient
H0     = simplify(subs(diff(diff(E,b_1),b_1),b_1,0),1000); % Hessian
E_quad = E0 + b_1*G0 + b_1^2*H0/2;                         % Local quadratic approximation

% Gradients (g) to use:
g = W(1,1)*vx_1 + mx_1*W(1,:)*(mx-mu) - 1
if simplify(G0 - g)~=0, disp('There''s a problem.'); end

% Hessian approximation (H) to use:
fprintf('g>0: ');
H  = W(1,1)*(mx_1^2+vx_1) + 1 + g
if simplify(H-H0)~=0, disp('There''s a problem.'); end
fprintf('g<0: ');
H  = W(1,1)*(mx_1^2+vx_1) + 1
Ha = simplify(subs(H0,mx_2,solve(G0==0,mx_2)),1000); % Check the workings
if simplify(H+g-H0)~=0, disp('There''s a problem.'); end

 *
 * nf   - Vector of dimensions (n_x, n_y, n_z, P).
 * mf   - E[f], dimensions nf.
 * vf   - Var[f], dimensions nf.
 * gmm  - Gaussian mixture model data structure.
 * nm   - Dimensions of log tissue priors (4 elements).
 * skip - Sampling density for log tissue priors (in x, y and z).
 * lkp  - Lookup table relating Gaussians to tissue classes.
 * lp   - Log tissue priors
 * g1   - Output gradients (n_x, n_y, n_z).
 * g2   - Output Hessian   (n_x, n_y, n_z).
 *
 */
static double INUgrads(mwSize nf[], float mf[], float vf[], unsigned char label[],
              mwSize K, GMMtype gmm[], double lnP[],
              mwSize nm[], mwSize skip[], mwSize lkp[], float lp[],
              mwSize index[],
              float  g1[], float g2[])
{
    mwSize P, Nf, Nm, i0,i1,i2, n0,n1,n2;
    double ll=0.0, mx[MaxChan], vx[MaxChan], p[128];

    P  = nf[3];
    Nf = nf[0]*nf[1]*nf[2];
    Nm = nm[0]*nm[1]*nm[2];

    n2 = nm[2]/skip[2]; if (n2>nf[2]) n2 = nf[2];
    n1 = nm[1]/skip[1]; if (n1>nf[1]) n1 = nf[1];
    n0 = nm[0]/skip[0]; if (n0>nf[0]) n0 = nf[0];

    if (P>=MaxChan || K>=128) return NAN;

    for(i2=0; i2<n2; i2++)
    {
        for(i1=0; i1<n1; i1++)
        {
            mwSize off_f, off_m;
            off_f = nf[0]*(i1         + nf[1]*i2);
            off_m = nm[0]*(i1*skip[1] + nm[1]*i2*skip[2]);
            for(i0=0; i0<n0; i0++)
            {
                mwSize i, im, code;
                i    = i0         + off_f;
                im   = i0*skip[0] + off_m;
                code = get_vox(Nf,P,mf+i,vf+i,mx,vx);
                if (code!=0 && logpriors(Nm, lp+im, K, lkp, p)!=0)
                {
                    if (label!=NULL) Dloglikelihoods((mwSize)(label[im]), K, lnP, p);
                    Nloglikelihoods(K, gmm, code, mx, vx, p);
                    ll += softmax1(K,p,p);
                    if (index[code]!=Undefined)
                    {
                        double g=0.0, h=0.0, *mu, *W, *nu;
                        mwSize Po = gmm[code].P, j, nc, k;
                        mu    = gmm[code].mu;
                        W     = gmm[code].W;
                        nu    = gmm[code].nu;
                        nc    = index[code];
                        for(k=0; k<K; k++)
                        {
                            double gk = 0.0, nup = nu[k]*p[k];
                            for(j=0; j<Po; j++)
                                gk += (mx[j]-mu[j+Po*k])*W[j+Po*(nc+Po*k)];
                            g += nup*gk;
                            h += nup*W[nc+Po*(nc+Po*k)];
                        }
                        g  = g*mx[nc]+h*vx[nc]        - 1.0;
                        h  = h*(mx[nc]*mx[nc]+vx[nc]) + 1.0;
                        h += fabs(g);

                        g1[i] = (float)g;
                        g2[i] = (float)h;
                    }
                }
            }
        }
    }
    return ll;
}

/* Construct a vector indicating which of the available
 * data correponds with the ic'th image
 *
 * P     - Number of imag volumes.
 * ic    - Index of image of interest.
 * index - Index of the ic't volume for each code.
 */
static void make_index(mwSize P, mwSize ic, mwSize index[])
{
    mwSize code,i,i1;
    for(code=0; code<(mwSize)1<<P; code++)
    {
        if ((code & (mwSize)1<<ic)!=0)
        {
            for(i=0,i1=0; i<ic; i++)
                if ((code & (mwSize)1<<i)!=0) i1++;
            index[code] = i1;
        }
        else
            index[code] = Undefined;
    }
}


/* Constructs a gmm structure from mu, b, W, nu, and gam, as well as a vector
 * of indices. These are then used when calling INUgrads.
 */
double call_INUgrads(mwSize nf[], float mf[], float vf[], unsigned char label[],
    mwSize K, double mu[], double b[], double W[], double nu[], double gam[], double lnP[],
    mwSize nm[], mwSize skip[], mwSize lkp[], float lp[],
    mwSize ic,
    float g1[], float g2[])
{
    mwSize P = nf[3];
    GMMtype *gmm;
    double ll;
    mwSize *index;

    if (P>=MaxChan || K>=128) return NAN;
    if ((gmm = sub_gmm(P, K, mu, b, W, nu, gam))==NULL) return NAN;
    index = (mwSize *)calloc((size_t)1<<(size_t)P, sizeof(mwSize));
    if (index == NULL)
    {
        (void)free((void *)gmm);
        return NAN;
    }
    make_index(P, ic, index);
    ll = INUgrads(nf, mf, vf, label, K, gmm, lnP, nm, skip, lkp, lp, index, g1, g2);
    (void)free((void *)gmm);
    (void)free((void *)index);
    return ll;
}

/* Fill in missing data, replacing NaNs with expectations
 *
 * nf   - Vector of dimensions (n_x, n_y, n_z, P).
 * mf   - E[f], dimensions nf.
 * vf   - Var[f], dimensions nf.
 * gmm  - Gaussian mixture model data structure.
 * nm   - Dimensions of log tissue priors (4 elements).
 * skip - Sampling density for log tissue priors (in x, y and z).
 * lkp  - Lookup table relating Gaussians to tissue classes.
 * lp   - Log tissue priors
 * mx1  - Output data (n_x, n_y, n_z, P).
 */
static int fill_missing(mwSize nf[], float mf[], float vf[], unsigned char label[],
              mwSize K, GMMtype gmm[], MissInfType missinf[],double lnP[],
              mwSize nm[], mwSize skip[], mwSize lkp[], float lp[],
              float  mx1[])
{
    mwSize P, Nf, Nm, i0,i1,i2, n0,n1,n2;
    double mx[MaxChan], vx[MaxChan], p[128];

    P  = nf[3];
    Nf = nf[0]*nf[1]*nf[2];
    Nm = nm[0]*nm[1]*nm[2];

    n2 = nm[2]/skip[2]; if (n2>nf[2]) n2 = nf[2];
    n1 = nm[1]/skip[1]; if (n1>nf[1]) n1 = nf[1];
    n0 = nm[0]/skip[0]; if (n0>nf[0]) n0 = nf[0];

    if (P>=MaxChan || K>=128) return -1;

    for(i2=0; i2<n2; i2++)
    {
        for(i1=0; i1<n1; i1++)
        {
            mwSize off_f, off_m;
            off_f = nf[0]*(i1         + nf[1]*i2);
            off_m = nm[0]*(i1*skip[1] + nm[1]*i2*skip[2]);
            for(i0=0; i0<n0; i0++)
            {
                mwSize i, im, code;
                i    = i0         + off_f;
                im   = i0*skip[0] + off_m;
                code = get_vox(Nf,P,mf+i,vf+i,mx,vx);
                if (logpriors(Nm, lp+im, K, lkp, p)!=0)
                {
                    double *mm, *mo, *Wt;
                    unsigned char *obs;
                    mwSize Po,Pm, j, jm,jo;

                    if (label!=NULL) Dloglikelihoods((mwSize)(label[i]), K, lnP, p);
                    Nloglikelihoods(K, gmm, code, mx, vx, p);
                    (void)softmax1(K,p,p);

                    Po    = missinf[code].Po;
                    Pm    = missinf[code].Pm;
                    mm    = missinf[code].m_m;
                    mo    = missinf[code].m_o;
                    Wt    = missinf[code].Wt;
                    obs   = missinf[code].obs;
                    for(j=0, jm=0, jo=0; j<P; j++)
                    {
                        if (obs[j]!=(unsigned char)0)
                        {   /* Simply copy existing data */
                            mx1[i+j*Nf] = mx[jo++];
                        }
                        else
                        {  /* Infer missing data using Bishop's PRML eq 2.96 & 2.97 */
                            mwSize k;
                            double tmp1;
                            for(k=0, tmp1=0.0; k<K; k++)
                            {
                                mwSize io;
                                double tmp = mm[jm+Pm*k];
                                for(io=0; io<Po; io++)
                                    tmp += Wt[jm+Pm*(io+Po*k)]*(mo[io+Po*k]-mx[io]);
                                tmp1 += tmp*p[k]; /* weight average by responsibilities */
                            }
                            mx1[i+j*Nf] = (float)tmp1;
                            jm ++;
                        }
                    }
                }
                else
                {
                    unsigned char *obs;
                    mwSize j, jo;
                    obs   = missinf[code].obs;
                    for(j=0, jo=0; j<P; j++)
                    {
                        if (obs[j]!=(unsigned char)0)
                            mx1[i+j*Nf] = mx[jo++];
                    }
                }
            }
        }
    }
    return 0;
}


int call_fill_missing(mwSize nf[], float mf[], float vf[], unsigned char label[],
    mwSize K, double mu[], double b[], double W[], double nu[], double gam[], double lnP[],
    mwSize nm[], mwSize skip[], mwSize lkp[], float lp[],
    float mx1[])
{
    mwSize P = nf[3];
    GMMtype *gmm;
    MissInfType *missinf;
    int sts;

    if (P>=MaxChan || K>=128) return -1;
    if ((gmm = sub_gmm(P, K, mu, b, W, nu, gam))==NULL) return -1;
    if ((missinf = prepare_missinf(P, K, W, mu))==NULL)
    {
        (void)free((void *)gmm);
        return -1;
    }

    sts = fill_missing(nf, mf, vf, label, K, gmm, missinf, lnP, nm, skip, lkp, lp, mx1);
    (void)free((void *)gmm);
    (void)free((void *)missinf);
    return sts;
}
