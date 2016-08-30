/*
 * $Id: spm_jsonread.c 6863 2016-08-30 14:56:27Z guillaume $
 * Guillaume Flandin
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "jsmn.h"
#include "mex.h"

static int create_struct(char *js, jsmntok_t *tok, mxArray **mx);

static char * get_string(char *js, int start, int end) {
    js[end] = '\0';
    return js + start;
}

static void valid_fieldname(char **field) {
    char *f = *field;
    while (f[0] != '\0') {
        if ( ((f[0] >= '0') && (f[0] <= '9')) 
          || ((f[0] >= 'a') && (f[0] <= 'z'))
          || ((f[0] >= 'A') && (f[0] <= 'Z'))) {
        }
        else {
        	f[0] = '_';
        }
        f++;
    }
    if ( (*field[0] == '_')
      || ((*field[0] >= '0') && (*field[0] <= '9')) ) {
        *field = *field - 1;
        *field[0] = 'x';
    }
}

static int primitive(char *js, jsmntok_t *tok, mxArray **mx) {
    mxArray *ma = NULL;
    int sts;
    switch (js[tok->start]) {
        case 't' :
            *mx =  mxCreateLogicalScalar(1);
            break;
        case 'f' :
            *mx =  mxCreateLogicalScalar(0);
            break;
        case 'n' :
            *mx =  mxCreateDoubleMatrix(0,0,mxREAL);
            break;
        default: /* '-', '0'..'9' */
            ma =  mxCreateString(get_string(js, tok->start, tok->end));
            sts = mexCallMATLAB(1, mx, 1, &ma, "str2double");
            if (sts != 0) {
                mexErrMsgTxt("Conversion from string to double failed.");
            }
            mxDestroyArray(ma);
            break;
    }
    return 1;
}

static int value(char *js, jsmntok_t *tok, mxArray **mx) {
    *mx = mxCreateString(get_string(js, tok->start, tok->end));
    return 1;
}

static int array(char *js, jsmntok_t *tok, mxArray **mx) {
    int i, j;
    mxArray *ma = NULL;
#ifndef MATLAB_MEX_FILE
    int sts;
#else
    mxArray *sts = NULL;
#endif
    mxArray *array[1];
    
    *mx = mxCreateCellMatrix(tok->size, 1);
    for (i = 0, j = 0; i < tok->size; i++) {
        j += create_struct(js, tok+1+j, &ma);
        mxSetCell(*mx, i, ma);
    }
    /* to do: nested arrays (of booleans, numeric, structure, strings)*/
    /* '[1,2]' => [1;2] */
    /* '[[1,2]]' => [1 2] */
    /* '[[1,2],[3,4]]' => [1 2;3 4] */
    
    /* Try to convert cell array into array*/
#ifndef MATLAB_MEX_FILE
    mexSetTrapFlag(1);
    sts = mexCallMATLAB(1, array, 1, mx, "cell2mat");
    mexPrintf("sts = %d\n",sts);
    if ((sts == 0) && (!mxIsChar(array[0]))) {
        return j+1; /* temporary always return cell arrays in Octave */
#else
    sts = mexCallMATLABWithTrap(1, array, 1, mx, "cell2mat");
    if ((sts == NULL) && (!mxIsChar(array[0]))) {
#endif
        mxDestroyArray(*mx);
        *mx = *array;
    }
    return j+1;
}

static int object(char *js, jsmntok_t *tok, mxArray **mx) {
    int i, j, k;
    mxArray *ma = NULL;
    char *field = NULL;
    if (tok->size == 0) {
        *mx = mxCreateStructMatrix(1, 1, 0, NULL);
        return 1;
    }
    for (i = 0, j = 0; i < tok->size; i++) {
        field = get_string(js, (tok+1+j)->start, (tok+1+j)->end); /* check it is a JSMN_STRING */
        valid_fieldname(&field);
        j++;
        if (i == 0) {
            *mx = mxCreateStructMatrix(1, 1, 1, (const char**)&field);
        }
        else {
            k = mxAddField(*mx, field);
            if (k == -1)
                mexErrMsgTxt("mxAddField()");
        }
        j += create_struct(js, tok+1+j, &ma);
        mxSetFieldByNumber(*mx, 0, i, ma);
    }
    return j+1;
}

static int create_struct(char *js, jsmntok_t *tok, mxArray **mx) {
    if (tok->type == JSMN_PRIMITIVE) {
        return primitive(js, tok, mx);
    } else if (tok->type == JSMN_STRING) {
        return value(js, tok, mx);
    } else if (tok->type == JSMN_OBJECT) {
        return object(js, tok, mx);
    } else if (tok->type == JSMN_ARRAY) {
        return array(js, tok, mx);
    }
    return 0;
}

static jsmntok_t * parse(const char *js, size_t jslen) {
    int r;
    jsmn_parser p;
    jsmntok_t *tok = NULL;
    size_t tokcount = 2;
    
    jsmn_init(&p);
    tok = mxMalloc(sizeof(*tok) * tokcount);
    if (tok == NULL) {
        mexErrMsgTxt("mxMalloc()");
    }
    
    for (;;) {
        r = jsmn_parse(&p, js, jslen, tok, tokcount);
        if (r < 0) {
            if (r == JSMN_ERROR_NOMEM) {
                tokcount = tokcount * 2;
                tok = mxRealloc(tok, sizeof(*tok) * tokcount);
                if (tok == NULL) {
                    mexErrMsgTxt("mxRealloc()");
                }
            }
            else if ((r == JSMN_ERROR_INVAL) || (r == JSMN_ERROR_PART)) {
                mexErrMsgTxt("Invalid or incomplete JSON.");
            }
            else {
                mexErrMsgTxt("Unknown JSON parsing error.");
            }
        }
        else {
            break;
        }
    }
    
    return tok;
}

static char * get_data(const mxArray * mx, size_t * jslen) {
    /* should attempt to minimise copy */
    int i, filename, sts;
    mxArray *ma = NULL;
    char *js = NULL;

    js = mxArrayToString(mx);
    if (js == NULL) {
        mexErrMsgTxt("mxArrayToString()");
    }
    *jslen = strlen(js);
    if (*jslen == 0)
        mexErrMsgTxt("Empty JSON.");
    
    /* detect whether input string is a filename */
    for (i = 0, filename = 1; i < *jslen; i++) {
        if ((js[i] == '{') || (js[i] == '[')) {
            filename = 0;
            break;
        }
    }
    if (filename == 1) {
        mxFree(js);
        sts = mexCallMATLAB(1, &ma, 1, (mxArray **)&mx, "fileread");
        if (sts != 0) {
            mexErrMsgTxt("Cannot read JSON file.");
        }
        js = mxArrayToString(ma);
        if (js == NULL) {
            mexErrMsgTxt("mxArrayToString()");
        }
        mxDestroyArray(ma);
        *jslen = strlen(js);
    }
    return js;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    char *js = NULL;
    size_t jslen = 0;
    jsmntok_t *tok = NULL;

    /* Validate input arguments */
    if (nrhs == 0) {
        mexErrMsgTxt("Not enough input arguments.");
    }
    else if (nrhs > 1) {
        mexErrMsgTxt("Too many input arguments.");
    }
    if (!mxIsChar(prhs[0])) {
        mexErrMsgTxt("Input must be a string.");
    }

    /* Get JSON data as char array */
    js = get_data(prhs[0], &jslen);

    /* Parse JSON data */
    tok = parse(js, jslen);

    /* Create output structure */
    create_struct(js, tok, &plhs[0]);

    mxFree(js);
    mxFree(tok);
}
