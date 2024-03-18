/*
 * Guillaume Flandin
 * Copyright (C) 2008-2023 Wellcome Centre for Human Neuroimaging
 */

#include <string.h>
#include "mex.h"

/*
 * https://stackoverflow.com/a/37109258
 * by polfosol: https://stackoverflow.com/users/5358284/polfosol
 */

static const char* B64chars = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

static const int B64index[256] =
{
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  62, 63, 62, 62, 63,
    52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 0,  0,  0,  0,  0,  0,
    0,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14,
    15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 0,  0,  0,  0,  63,
    0,  26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
    41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51
};

unsigned char * b64encode(const unsigned char*p, size_t *len)
{
    size_t l = (*len + 2) / 3 * 4;
    unsigned char *result = mxMalloc(l); /* '=' */
    size_t i, j = 0, pad = *len % 3;
    const size_t last = *len - pad;
    int n;
    result[l - 1] = '=';
    result[l - 2] = '=';

    for (i = 0; i < last; i += 3)
    {
        n = (int)(p[i]) << 16 | (int)(p[i + 1]) << 8 | p[i + 2];
        result[j++] = B64chars[n >> 18];
        result[j++] = B64chars[n >> 12 & 0x3F];
        result[j++] = B64chars[n >> 6 & 0x3F];
        result[j++] = B64chars[n & 0x3F];
    }
    if (pad) 
    {
        n = --pad ? (int)(p[last]) << 8 | p[last + 1] : p[last];
        result[j++] = B64chars[pad ? n >> 10 & 0x3F : n >> 2];
        result[j++] = B64chars[pad ? n >> 4 & 0x03F : n << 4 & 0x3F];
        result[j++] = pad ? B64chars[n << 2 & 0x3F] : '=';
    }
    *len = l;
    return result;
}

unsigned char * b64decode(const unsigned char *p, size_t *len)
{
    int n;
    size_t i, j = 0,
        pad1 = *len % 4 || p[*len - 1] == '=',
        pad2 = pad1 && (*len % 4 > 2 || p[*len - 2] != '=');
    const size_t last = (*len - pad1) / 4 << 2;
    unsigned char *result = mxMalloc(last / 4 * 3 + pad1 + pad2); /* '\0' */
    unsigned char *str = (unsigned char*) &result[0];

    if (len == 0) return NULL;
    
    for (i = 0; i < last; i += 4)
    {
        n = B64index[p[i]]     << 18 | 
            B64index[p[i + 1]] << 12 | 
            B64index[p[i + 2]] << 6  | 
            B64index[p[i + 3]];
        str[j++] = n >> 16;
        str[j++] = n >> 8 & 0xFF;
        str[j++] = n & 0xFF;
    }
    if (pad1)
    {
        n = B64index[p[last]] << 18 | B64index[p[last + 1]] << 12;
        str[j++] = n >> 16;
        if (pad2)
        {
            n |= B64index[p[last + 2]] << 6;
            str[j++] = n >> 8 & 0xFF;
        }
    }
    *len = last / 4 * 3 + pad1 + pad2;
    return result;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    char *action = NULL;
    unsigned char *str = NULL;
    size_t n;
    
    /* Check for proper number of arguments */
    if (nrhs < 2)
        mexErrMsgTxt("Not enough input arguments.");
    else if (nrhs > 2)
        mexErrMsgTxt("Too many input arguments.");
    else if (nlhs > 1)
        mexErrMsgTxt("Too many output arguments.");
    
    /* The first input must be a string */
    if (!mxIsChar(prhs[0]))
        mexErrMsgTxt("First input must be a string.");
    action = mxArrayToString(prhs[0]);
    
    /* The second input must be a real uint8 array */
    if (!mxIsUint8(prhs[1]) || mxIsComplex(prhs[1]))
        mexErrMsgTxt("Second input must be a real uint8 array.");
    
    /* Base64 encoding/decoding */
    n = mxGetNumberOfElements(prhs[1]);
    if (!strcmp(action,"encode")) {
        str = b64encode((unsigned char*)mxGetPr(prhs[1]),&n);
    }
    else if (!strcmp(action,"decode")) {
        str = b64decode((unsigned char*)mxGetPr(prhs[1]),&n);
    }
    else
        mexErrMsgTxt("Unknown action.");
    
    /* Return output as uint8 row vector */
    /*plhs[0] = mxCreateNumericMatrix(1,n,mxUINT8_CLASS,mxREAL);
    memcpy(mxGetData(plhs[0]), str, n);*/
    plhs[0] = mxCreateNumericMatrix(0,0,mxUINT8_CLASS,mxREAL);
    mxFree(mxGetPr(plhs[0]));
    mxSetPr(plhs[0],(double *)str);
    mxSetM(plhs[0],1);
    mxSetN(plhs[0],n);
    
    mxFree(action);
}
