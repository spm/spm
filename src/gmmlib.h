/* $Id$ */
/* (c) John Ashburner, Mikael Brudfors & Yael Balbastre (2020) */
#include<math.h>
#include<stdlib.h>

void choldc(unsigned long long n, double a[], double p[]);
void cholls(unsigned long long n, const double a[], const double p[], const double b[], double x[]);
void test(size_t M, size_t K, double *mu, double *b, double *W, double *nu, double *gam);
double psi(double z);

void space_needed(size_t M, size_t K, size_t *m0, size_t *m1, size_t *m2);

double call_suffstats_missing(size_t nf[], float mf[], float vf[],
    size_t K, double mu[], double b[], double W[], double nu[], double gam[],
    size_t nm[], size_t skip[], size_t lkp[], float lp[],
    double s0_ptr[], double s1_ptr[], double s2_ptr[]);

double call_responsibilities(size_t nf[], size_t skip[], float mf[], float vf[],
    size_t K, double mu[], double b[], double W[], double nu[], double gam[],
    size_t K1, size_t lkp[], float lp[],
    float r[]);

double call_INUgrads(size_t nf[], float mf[], float vf[],
    size_t K, double mu[], double b[], double W[], double nu[], double gam[],
    size_t nm[], size_t skip[], size_t lkp[], float lp[],
    size_t ic, float fc[],
    float g1[], float g2[]);

