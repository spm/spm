/*
 * John Ashburner & Matthew Brett
 * Copyright (C) 1999-2022 Wellcome Centre for Human Neuroimaging
 */

/* Routines for accessing datatypes for images */

#ifndef _SPM_GETDATA_H_
#define _SPM_GETDATA_H_

signed   short     getshort( signed   short x);
unsigned short     getushort(unsigned short x);
signed   int       getint(   signed   int x);
unsigned int       getuint(  unsigned int x);
signed   long long getint64( signed   long long x);
unsigned long long getuint64(unsigned long long x);
float  getfloat( float  x);
double getdouble(double x);

#endif /* _SPM_GETDATA_H_ */
