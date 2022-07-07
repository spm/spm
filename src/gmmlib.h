/*
 * John Ashburner, Mikael Brudfors & Yael Balbastre
 * Copyright (C) 2020-2022 Wellcome Centre for Human Neuroimaging
 */

#include<math.h>
#include<stdlib.h>

void choldc(mwSize n, double a[], double p[]);
void cholls(mwSize n, const double a[], const double p[], const double b[], double x[]);
void test(mwSize M, mwSize K, double *mu, double *b, double *W, double *nu, double *gam);
double psi(double z);

void space_needed(mwSize M, mwSize K, mwSize *m0, mwSize *m1, mwSize *m2);

double call_suffstats_missing(mwSize nf[], float mf[], float vf[], unsigned char label[],
    mwSize K, double mu[], double b[], double W[], double nu[], double gam[], double lnP[],
    mwSize nm[], mwSize skip[], mwSize lkp[], float lp[],
    double s0_ptr[], double s1_ptr[], double s2_ptr[], double H[]);

double call_responsibilities(mwSize nf[], mwSize skip[], float mf[], float vf[], unsigned char label[],
    mwSize K, double mu[], double b[], double W[], double nu[], double gam[], double lnP[],
    mwSize K1, mwSize lkp[], float lp[],
    float r[]);

double call_INUgrads(mwSize nf[], float mf[], float vf[], unsigned char label[],
    mwSize K, double mu[], double b[], double W[], double nu[], double gam[], double lnP[],
    mwSize nm[], mwSize skip[], mwSize lkp[], float lp[],
    mwSize ic,
    float g1[], float g2[]);

int    call_fill_missing(mwSize nf[], float mf[], float vf[], unsigned char label[],
    mwSize K, double mu[], double b[], double W[], double nu[], double gam[], double lnP[],
    mwSize nm[], mwSize skip[], mwSize lkp[], float lp[],
    float mf1[]);
