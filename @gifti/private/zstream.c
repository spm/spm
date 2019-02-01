/*
 * $Id: zstream.c 7523 2019-02-01 11:31:08Z guillaume $
 * Guillaume Flandin
 */

/* mex -O CFLAGS='$CFLAGS -std=c99' -largeArrayDims zstream.c */

/* setenv CFLAGS "`mkoctfile -p CFLAGS` -std=c99" */
/* mkoctfile --mex zstream.c */

/* miniz: https://github.com/richgel999/miniz */
#define MINIZ_NO_STDIO
#define MINIZ_NO_ARCHIVE_APIS
#define MINIZ_NO_TIME
#define MINIZ_NO_ZLIB_APIS
#include "miniz.c"

#include "mex.h"

/* --- GATEWAY FUNCTION --- */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

char *action = NULL;
unsigned char *IN = NULL, *OUT = NULL;
size_t INlen, OUTlen;
int flag = 0;

/* Check for proper number of arguments */
if (nrhs < 2)
    mexErrMsgTxt("Not enough input arguments.");
else if (nrhs > 2)
    mexErrMsgTxt("Too many input arguments.");
else if (nlhs > 1)
    mexErrMsgTxt("Too many output arguments.");

/* The input ACTION must be a string */
if (!mxIsChar(prhs[0]))
    mexErrMsgTxt("Input ACTION must be a string.");
action = mxArrayToString(prhs[0]);

/* The input IN must be a real uint8 array */
if (!mxIsUint8(prhs[1]) || mxIsComplex(prhs[1]))
    mexErrMsgTxt("Input IN must be a real uint8 array.");

INlen = mxGetNumberOfElements(prhs[1]);
IN = mxGetData(prhs[1]);

/* zlib stream (zlib header with adler32 checksum) or raw deflate stream */
if (!strcmp(action,"D")) {
    flag = TINFL_FLAG_PARSE_ZLIB_HEADER;
}
else if (!strcmp(action,"C")) {
    flag = TDEFL_WRITE_ZLIB_HEADER;
}
    
if (!strcmp(action,"D") || !strcmp(action,"d")) {

    /* Decompress data */
    OUT = tinfl_decompress_mem_to_heap(IN, INlen, &OUTlen, flag);

    if (OUT == NULL)
        mexErrMsgTxt("Error when decompressing data.");
}
else if (!strcmp(action,"C") || !strcmp(action,"c")) {
    /* Compress data */
    OUT = tdefl_compress_mem_to_heap(IN, INlen, &OUTlen, flag);
    
    if (OUT == NULL)
        mexErrMsgTxt("Error when compressing data.");
}
else {
    mexErrMsgTxt("Unknown ACTION type.");
}

/* Store output */
plhs[0] = mxCreateNumericMatrix(OUTlen,1,mxUINT8_CLASS,mxREAL);
if (plhs[0] == NULL)
    mexErrMsgTxt("Error when creating output variable.");

memcpy(mxGetData(plhs[0]), OUT, OUTlen);

mxFree(action);
mz_free(OUT);

}
