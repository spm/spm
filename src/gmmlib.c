/*
 * $Id: gmmlib.c 8014 2020-11-23 18:25:11Z john $
 * John Ashburner, Mikael Brudfors & Yael Balbastre
 */

#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#define EXP(x) fastexp(x)
typedef signed long long int  Int64;

typedef struct
{
    double *s2;
    double *s1;
    double *s0;
} SStype;

typedef struct
{
    size_t  P;
    double *mu;
    double *b;
    double *W;
    double *nu;
    double *gam;
    double *conN;
    double *conT;
} GMMtype;

static const double pi = 3.1415926535897931;
static const size_t MaxChan=(size_t)50; /* largest integer valued float is 2^52 */
static const size_t Undefined=(size_t)0xFFFFFFFFFFFFF;

static double fastexp(double x)
{
    double r, rr;
    signed long long i;
    static double lkp_mem[256], *exp_lkp = lkp_mem+128;

    /* exp(i+r) = exp(i)*exp(r), where:
     *     exp(i) is from the lookup table;
     *     exp(r) is from a generalised continued fraction
     *            https://en.wikipedia.org/wiki/Exponential_function#Continued_fractions_for_ex
     *
     * Should not encounter values more extreme than -128 or 127,
     * particularly as the upper limit of x will be 0 and values
     * of x below log(eps)=-36.04 should be numerically equivalent.*/
    i  = (signed long long)rint(x);
    if (i<-128) i = -128;
    if (i> 127) i =  127;
    if (exp_lkp[i]==0.0) exp_lkp[i] = exp((double)i);

    r  = x - (double)i;
    rr = r*r;
/*  return exp_lkp[i] * (1.0+2.0*r/(2.0-r+rr/(6.0+rr/(10.0+rr/14.0))));
 *  return exp_lkp[i] * (1.0+2.0*r/(2.0-r+rr/(6.0+rr/(10.0)))); */
    return exp_lkp[i] * (1.0+2.0*r/(2.0-r+rr/6.0));
}


SStype *suffstat_pointers(size_t P, size_t K, double *s0_ptr, double *s1_ptr, double *s2_ptr)
{
    SStype /*@NULL@*/ *suffstat;
    suffstat = (SStype *)calloc((size_t)1<<P,sizeof(SStype));
    if (suffstat != NULL)
    {
        size_t code,o0=0,o1=0,o2=0;
        for(code=0; code<((size_t)1<<P); code++)
        {
            size_t i, Po;
            for(i=0, Po=0; i<P; i++) Po += (code>>i) & (size_t)1;
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


static size_t get_vox(size_t N1, size_t P, float mf[], float vf[], /*@out@*/ double x[], /*@out@*/ double v[])
{
    size_t j, j1, o, code;
    for(j=0, j1=0, o=0, code=0; j<P; j++, o+=N1)
    {
        double tmp = (double)mf[o];
        if (isfinite(tmp))
        {
            x[j1] = tmp;
            v[j1] = vf[o];
            code |= (size_t)1<<j;
            j1++;
        }
    }
    return code;
}

static double lse(size_t K, double q[])
{
    size_t k;
    double mx, s;
    for(k=1, mx=q[0]; k<K; k++)
    {
        if (q[k]>mx) mx = q[k];
    }
    for(k=0, s=0.0; k<K; k++)
        s += EXP(q[k]-mx);
    return log(s)+mx;
}


static double softmax1(size_t K, double q[], /*@out@*/ double p[])
{
    size_t k;
    double mx, s;
    for(k=1, mx=q[0]; k<K; k++)
        if (q[k]>mx) mx = q[k];
    for(k=0, s=0.0; k<K; k++)
        s += (p[k] = EXP(q[k]-mx));
    for(k=0; k<K; k++)
        p[k] /= s;
    return log(s) + mx;
}


static double softmax(size_t K, double q[], /*@out@*/ double p[])
{
    size_t k;
    double mx, s;
    for(k=0, mx=0; k<K; k++)
        if (q[k]>mx) mx = q[k];
    for(k=0, s=EXP(-mx); k<K; k++)
        s += (p[k] = EXP(q[k]-mx));
    for(k=0; k<K; k++)
        p[k] /= s;
    return log(s) + mx;
}


static double del2(size_t P, double mu[], double W[], double x[], double v[])
{
    size_t j,i;
    double d=0.0, r, *wj;
    for(j=0,wj=W; j<P; j++, wj+=P)
    {
        r  = x[j]-mu[j];
        d += wj[j]*(r*r+v[j]);
        for(i=j+1; i<P; i++)
            d += 2.0*wj[i]*r*(x[i]-mu[i]);
    }
    return d;
}


static double psi(double z)
{
    /* From http://web.science.mq.edu.au/~mjohnson/code/digamma.c */
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


static double Nresp(size_t K, GMMtype gmm[], size_t code, double x[], double v[], double p[])
{
    size_t P, k;
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
    return softmax1(K,p,p);
}


static double Tresp(size_t K, GMMtype gmm[], size_t code, double x[], double v[], double p[])
{
    size_t P, k;
    double *mu, *b, *W, *nu, *con;
    /*
       Compute other responsibilities from a mixture of Student's t distributions.
       See Eqns. 10.78-10.82 & B.68-B.72 in Bishop's PRML book.
       In practice, it only improves probabilities by a tiny amount.

       ln St(x|mu,Lam,tau)
       lgamma((tau + P)/2.0) - lgamma(tau/2.0) + sum(log(diag(chol(Lam)))) - (P/2.0)*log(tau*pi) -
       ((tau+P)/2)*log(1 + ((x-mu)'*Lam*(x-mu) + sum(diag(Lam).*vf))/tau)
       where:
           Lam = (nu+1-P)*beta/(1+beta)*W
           tau = nu+1-P

       lgamma((nu+1)/2.0) - lgamma((nu+1-P)/2.0) +  sum(log(diag(chol(Lam)))) -(P/2.0)*log((nu+1-P)*pi) -
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
    return softmax1(K,p,p);
}

static int get_priors(size_t N1, float *lp, size_t K, size_t *lkp, double *p)
{
    size_t k;
    double l;
    for(k=0; k<K; k++)
    {
        double lpk;
        lpk = (double)lp[N1*lkp[k]];
        if (!isfinite(lpk))
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


static void choldc(size_t n, double a[], /*@out@*/ double p[])
{
    Int64 i, j, k;
    double sm, sm0;

    sm0  = 1e-40;
    for(i=0; i<(Int64)n; i++) sm0 = sm0 + a[i*n+i];
    sm0 *= 1e-7;
    sm0 *= sm0;

    for(i=0; i<(Int64)n; i++)
    {
        for(j=i; j<(Int64)n; j++)
        {
            sm = a[i*n+j];
            for(k=i-1; k>=0; k--)
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


static void cholls(size_t n, const double a[], const double p[],
            const double b[], /*@out@*/ double x[])
{
    Int64 i, k;
    double sm;

    for(i=0; i<(Int64)n; i++)
    {
        sm = b[i];
        for(k=i-1; k>=0; k--)
            sm -= a[i*n+k]*x[k];
        x[i] = sm/p[i];
    }
    for(i=(Int64)n-1; i>=0; i--)
    {
        sm = x[i];
        for(k=i+1; k<(Int64)n; k++)
            sm -= a[k*n+i]*x[k];
        x[i] = sm/p[i];
    }
}


static size_t factorial(size_t n)
{
    static size_t products[21];
    if (products[0]==0)
    {
        size_t i;
        products[0] = 1;
        for(i=1; i<21; i++)
            products[i] = products[i-1]*i;
    }
    return products[n];
}

void space_needed(size_t P, size_t K, size_t *m0, size_t *m1, size_t *m2)
{
    size_t m;
    for(m=0, *m0=0, *m1=0, *m2=0; m<=P; m++)
    {
        size_t nel;
        nel = K*factorial(P)/(factorial(m)*factorial(P - m));
        *m0 += nel;
        *m1 += nel*m;
        *m2 += nel*m*m;
    }
}

static GMMtype *allocate_gmm(size_t P, size_t K)
{
    size_t o, code, i, n0=0,n1=0,n2=0;
    double *buf;
    unsigned char *bytes;
    GMMtype /*@NULL@*/ *gmm;
    space_needed(P, K, &n0, &n1, &n2);

    o     = ((size_t)1<<P)*sizeof(GMMtype);
    bytes = calloc(o+(n0*(size_t)5+n1+n2)*sizeof(double),1);
    gmm   = (GMMtype *)bytes;
    if (gmm!=NULL)
    {
        buf   = (double *)(bytes + o);
        o     = 0;
        for(code=0; code<((size_t)1<<P); code++)
        {
            size_t nel = 0;
            for(i=0; i<code; i++) nel += (code>>i) & 1;
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


static double invert(size_t P, double *W /* P*P */, double *S /* P*P */, double *T /* P*(P+1) */)
{
    size_t i, j, PP=P*P;
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


static GMMtype *sub_gmm(size_t P, size_t K, double *mu, double *b, double *W, double *nu, double *gam)
{
    const double log2pi = log(2*pi), log2 = log(2.0);
    double *S, *Si;
    GMMtype *gmm;
    size_t k, code, PP = P*P;

    if ((gmm = allocate_gmm(P,K)) == NULL) return gmm;
    if ((S   = (double *)calloc(P*((size_t)3*P+(size_t)1),sizeof(double))) == NULL)
    {
       (void)free((void *)gmm);
       return NULL;
    }
    Si = S + PP;
    for(k=0; k<K; k++)
    {
        double lgam = log(gam[k]);
        (void)invert(P,W+PP*k,S,S+PP);
        for(code=0; code<(size_t)1<<P; code++)
        {
            size_t j,j1, Po;
            double ld, eld;
            Po               = gmm[code].P;
            gmm[code].nu[k]  = nu[k] - (P-Po);
            gmm[code].b[k]   = b[k];
            gmm[code].gam[k] = lgam;
            for(j=0, j1=0; j<P; j++)
            {
                if ((((size_t)1<<j) & code) != 0)
                {
                    size_t i, i1;
                    gmm[code].mu[j1+Po*k] = mu[j+P*k];
                    for(i=0, i1=0; i<P; i++)
                    {
                        if ((((size_t)1<<i) & code) != 0)
                        {
                            Si[i1+Po*j1] = S[i+P*j];
                            i1++;
                        }
                    }
                    j1++;
                }
            }
            ld = invert(Po,Si,gmm[code].W+k*Po*Po,Si+Po*Po);

            /* Constant term for mixture of Gaussians
               E[ln N(x | m, L^{-1})] w.r.t. Gaussian-Wishart */
            for(j=0,eld=0.0; j<Po; j++) eld += psi((gmm[code].nu[k]-(double)j)*0.5);
            eld    += Po*log2 + ld;
            gmm[code].conN[k] = 0.5*(eld - Po*(log2pi+1.0/gmm[code].b[k])) + lgam;

            /* Constant term for mixture of T distributions */
            gmm[code].conT[k] = lgamma(0.5*(gmm[code].nu[k]+1.0)) - lgamma(0.5*(gmm[code].nu[k]+1.0-Po)) +
                                0.5*ld - 0.5*P*log((gmm[code].nu[k]+1-Po)*pi) + lgam;
        }
    }
    (void)free((void *)S);
    return gmm;
}


static double suffstats_missing(size_t nf[], float mf[], float vf[],
                      size_t K, GMMtype gmm[],
                      size_t nm[], size_t skip[], size_t lkp[], float lp[],
                      SStype suffstat[])
{
    size_t K1, i0,i1,i2, n2,n1,n0, P, Nf, Nm, code;
    double ll = 0.0, x[MaxChan], e[MaxChan], p[128];

    P  = nf[3];
    Nf = nf[0]*nf[1]*nf[2];
    K1 = nm[3];
    Nm = nm[0]*nm[1]*nm[2];

    n2 = nm[2]/skip[2]; if (n2>nf[2]) n2 = nf[2];
    n1 = nm[1]/skip[1]; if (n1>nf[1]) n1 = nf[1];
    n0 = nm[0]/skip[0]; if (n0>nf[0]) n0 = nf[0];

    for(i2=0; i2<n2; i2+=skip[2])
    {
        for(i1=0; i1<n1; i1+=skip[1])
        {
            size_t off_f, off_m;
            off_f = nf[0]*(i1         + nf[1]*i2);
            off_m = nm[0]*(i1*skip[1] + nm[1]*i2*skip[2]);
            for(i0=0; i0<n0; i0+=skip[0])
            {
                size_t i, im;
                i    = i0         + off_f;
                im   = i0*skip[0] + off_m;
                code = get_vox(Nf,P,mf+i,vf+i,x,e);
                if (code>0 && get_priors(Nm, lp+im, K, lkp, p)!=0)
                {
                    size_t j, j1, k, Po;
                    double *s0, *s1, *s2;
                    ll += Nresp(K, gmm, code, x, e, p);
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
                            double xj = x[j], px = pk*xj;
                            s1[j]      += px;
                            s2[j+Po*j] += pk*(xj*xj+e[j]);
                            for(j1=j+1; j1<Po; j1++)
                                s2[j1+Po*j] += px*x[j1];
                        }
                    }
                }
            }
        }
    }

    /* Add in upper triangle second order sufficiant statistics */
    for(code=1; code<((size_t)1<<P); code++)
    {
        size_t j, j1, k, Po;
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


double call_suffstats_missing(size_t nf[], float mf[], float vf[],
    size_t K, double mu[], double b[], double W[], double nu[], double gam[],
    size_t nm[], size_t skip[], size_t lkp[], float lp[],
    double s0_ptr[], double s1_ptr[], double s2_ptr[])
{
    size_t P = nf[3];
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
    ll = suffstats_missing(nf, mf, vf, K, gmm, nm, skip, lkp, lp, suffstat);
    (void)free((void *)gmm);
    (void)free((void *)suffstat);
    return ll;
}


static double responsibilities(size_t nf[], size_t skip[], float mf[], float vf[],
              size_t K, GMMtype *gmm,
              size_t K1, size_t lkp[], float lp[],
              float r[])
{
    size_t P, N1, i0,i1,i2;
    double ll = 0.0, x[MaxChan], e[MaxChan], p[128];

    P  = nf[3];
    N1 = nf[0]*nf[1]*nf[2];

    for(i2=0; i2<nf[2]; i2++)
    {
        for(i1=0; i1<nf[1]; i1++)
        {
            size_t off_f;
            off_f = nf[0]*(i1 + nf[1]*i2);
            for(i0=0; i0<nf[0]; i0++)
            {
                size_t i, code, k, k1;
                i    = i0+off_f;
                code = get_vox(N1,P,mf+i,vf+i,x,e);
                if (get_priors(N1, lp+i, K, lkp, p)!=0)
                {
                    if (code!=0)
                    {
                        if ((i2%skip[2])==0 && ((i1%skip[1])==0) & ((i0%skip[0])==0))
                            ll += Nresp(K, gmm, code, x, e, p);
                        else
                            ll += Tresp(K, gmm, code, x, e, p);
                        for(k=0; k<K; k++)
                        {
                            k1 = lkp[k];
                            if (k1<K1-1)
                                r[i+k1*N1] += p[k];
                        }
                    }
                    else
                    {
                        (void)softmax(K1,p,p);
                        for(k1=0; k1<K1-1; k1++)
                            r[i+k1*N1]  = NAN;
                        /*  r[i+k1*N1] += p[k1]; */
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

double call_responsibilities(size_t nf[], size_t skip[], float mf[], float vf[],
    size_t K, double mu[], double b[], double W[], double nu[], double gam[],
    size_t K1, size_t lkp[], float lp[],
    float r[])
{
    size_t P = nf[3];
    GMMtype *gmm;
    double ll;

    if (P>=MaxChan || K>=128) return NAN;
    if ((gmm      = sub_gmm(P, K, mu, b, W, nu, gam))==NULL) return NAN;

    ll = responsibilities(nf, skip, mf, vf, K, gmm, K1, lkp, lp, r);

    (void)free((void *)gmm);
    return ll;
}


static double INUgrads(size_t nf[], float mf[], float vf[],
              size_t K, GMMtype gmm[],
              size_t nm[], size_t skip[], size_t lkp[], float lp[],
              size_t index[], float fc[],
              float  g1[], float g2[])
{
    size_t P, Nf, Nm, i0,i1,i2, n0,n1,n2;
    double ll=0.0, x[MaxChan], e[MaxChan], p[128];

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
            size_t off_f, off_m;
            off_f = nf[0]*(i1         + nf[1]*i2);
            off_m = nm[0]*(i1*skip[1] + nm[1]*i2*skip[2]);
            for(i0=0; i0<n0; i0++)
            {
                size_t i, im, code;
                i    = i0         + off_f;
                im   = i0*skip[0] + off_m;
                code = get_vox(Nf,P,mf+i,vf+i,x,e);
                if (code!=0 && get_priors(Nm, lp+im, K, lkp, p)!=0)
                {
                    ll += Nresp(K, gmm, code, x, e, p);
                    if (index[code]!=Undefined)
                    {
                        double g=0.0, h=0.0, *mu, *W, *nu;
                        size_t Po = gmm[code].P, j, nc, k;
                        mu    = gmm[code].mu;
                        W     = gmm[code].W;
                        nu    = gmm[code].nu;
                        nc    = index[code];
                        x[nc] = fc[i];       /* replace E[f] with the point estimate */
                        for(k=0; k<K; k++)
                        {
                            double gk = 0.0, nup = nu[k]*p[k];
                            for(j=0; j<Po; j++)
                                gk += (x[j]-mu[j+Po*k])*W[j+Po*(nc+Po*k)];
                            g += nup*gk;
                            h += nup*W[nc+Po*(nc+Po*k)];
                        }
                        g = g*fc[i]       - 1.0;
                        h = h*fc[i]*fc[i] + 1.0;
                        if (g>0.0) h += g;

                        g1[i] = (float)g;
                        g2[i] = (float)h;
                    }
                }
            }
        }
    }
    return ll;
}


static void make_index(size_t P, size_t ic, size_t index[])
{
    size_t code,i,i1;
    for(code=0; code<(size_t)1<<P; code++)
    {
        if ((code & (size_t)1<<ic)!=0)
        {
            for(i=0,i1=0; i<ic; i++)
                if ((code & (size_t)1<<i)!=0) i1++;
            index[code] = i1;
        }
        else
            index[code] = Undefined;
    }
}

double call_INUgrads(size_t nf[], float mf[], float vf[],
    size_t K, double mu[], double b[], double W[], double nu[], double gam[],
    size_t nm[], size_t skip[], size_t lkp[], float lp[],
    size_t ic, float fc[],
    float g1[], float g2[])
{
    size_t P = nf[3];
    GMMtype *gmm;
    double ll;
    size_t *index;

    if (P>=MaxChan || K>=128) return NAN;
    if ((gmm = sub_gmm(P, K, mu, b, W, nu, gam))==NULL) return NAN;
    index = (size_t *)calloc((size_t)1<<P, sizeof(size_t));
    if (index == NULL)
    {
        (void)free((void *)gmm);
        return NAN;
    }
    make_index(P, ic, index);
    ll = INUgrads(nf, mf, vf, K, gmm, nm, skip, lkp, lp, index, fc, g1, g2);
    (void)free((void *)gmm);
    (void)free((void *)index);
    return ll;
}

