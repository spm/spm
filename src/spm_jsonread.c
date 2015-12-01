/*
 * $Id: spm_jsonread.c 6618 2015-12-01 16:25:38Z spm $
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
            *mx =  mxCreateDoubleScalar(mxGetNaN());
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
    *mx = mxCreateCellMatrix(tok->size, 1);
    for (i = 0; i < tok->size; i++) {
        j += create_struct(js, tok+1+j, &ma);
        mxSetCell(*mx, i, ma);
    }
    /* Call cell2mat if all cell elements are primitives (mxIsNumeric(ma))? */
    return j+1;
}

static int object(char *js, jsmntok_t *tok, mxArray **mx) {
    int i, j, k;
    mxArray *ma = NULL;
    char *field = NULL;
    for (i = 0, j = 0; i < tok->size; i++) {
        field = get_string(js, (tok+1+j)->start, (tok+1+j)->end); /* check it is a JSMN_STRING */
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
