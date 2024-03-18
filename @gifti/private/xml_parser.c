/*
 * Guillaume Flandin
 * Copyright (C) 2008-2023 Wellcome Centre for Human Neuroimaging
 */

#ifndef MATLAB_MEX_FILE
# undef  _LARGEFILE64_SOURCE
# define _LARGEFILE64_SOURCE
# include <stdio.h>
# include <sys/stat.h>
# if defined(__APPLE__)
#  define structStat struct stat
#  define getFileFstat fstat
# elif defined(SPM_WIN32)
#  define structStat struct _stati64
#  define getFileFstat _fstati64
# else
#  define structStat struct stat64
#  define getFileFstat fstat64
# endif
#else
# include "io64.h"
#endif
#include "mex.h"
#include "yxml.h"
#include <string.h>

/* mex CFLAGS='$CFLAGS -std=c99' xml_parser.c yxml.c */

static char * get_data(const mxArray * mx) {
  char *xml = NULL;
  int filename;
  size_t i, xmllen;

  xml = mxArrayToString(mx);
  if (xml == NULL)
    mexErrMsgTxt("mxArrayToString() failed");
  xmllen = strlen(xml);
  if (xmllen == 0)
    mexErrMsgTxt("Empty XML document");

  /* detect whether input string is a filename */
  for (i = 0, filename = 1; i < xmllen; i++)
    if (xml[i] == '<') {
      filename = 0;
      break;
    }
  if (filename == 1) {
    structStat statbuf;
    int64_T fileSize = 0;
    FILE *fp = fopen(xml, "rb");
    if (fp == NULL)
      mexErrMsgTxt("Cannot read XML document");
    if (getFileFstat(fileno(fp), &statbuf) == 0)
      fileSize = statbuf.st_size;
    else
      mexErrMsgTxt("Cannot get size of XML document");
    mxFree(xml);
    xml = mxMalloc(fileSize + 1);
    if (xml) {
      xml[fileSize] = '\0';
      if (fread (xml, 1, fileSize, fp) != (size_t)fileSize)
        mexErrMsgTxt("Cannot read the entire XML document");
    }
    else
      mexErrMsgTxt("Cannot allocate memory for the XML document");
    fclose(fp);
  }
  return xml;
}

static void resize_struct(mxArray * mx, mwSize N) {
  mwIndex i;
  int j, nfields = mxGetNumberOfFields(mx);
  size_t n = mxGetNumberOfElements(mx);
  mxSetData(mx, mxRealloc(mxGetData(mx), N * nfields * sizeof(mxArray *)));
  mxSetN(mx, N);
  for(i = n; i < N; i++)
    for(j = 0; j < nfields; j++)
      mxSetFieldByNumber(mx, i, j, NULL);
}

static inline void concatenate_array(mxArray * mx, int c, int n, int val) {
  mxArray *a = mxGetFieldByNumber(mx,c,n);
  if (a == NULL)
    mxSetFieldByNumber(mx, c, n, mxCreateDoubleScalar(val));
  else {
    size_t noe = mxGetNumberOfElements(a);
    mxSetPr(a, mxRealloc(mxGetPr(a), (noe + 1) * sizeof(double)));
    mxSetN(a, noe + 1);
    mxGetPr(a)[noe] = val;
    mxSetFieldByNumber(mx, c, n, a);
  }
}

static char XML_STACK[8*1024];

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  yxml_t x[1];
  char *xml = NULL, *doc = NULL;
  size_t nattr = 0;
  mxArray *mx = NULL;
  int XML_UID = 0;
  int XML_CURR_NODE = 0;
  char *XML_START = NULL;
  char *XML_CURR = NULL;
  /*
  const char **nodefields = (const char *[]){"type", "value", "attributes", "children", "uid", "parent"};
  const char **attrfields = (const char *[]){"key", "value"};
   */
  const char **nodefields = mxMalloc(sizeof(char*) * 6);
  const char **attrfields = mxMalloc(sizeof(char*) * 2);
  nodefields[0] = "type";
  nodefields[1] = "value";
  nodefields[2] = "attributes";
  nodefields[3] = "children";
  nodefields[4] = "uid";
  nodefields[5] = "parent";
  attrfields[0] = "key";
  attrfields[1] = "value";
  
  /* Initialisation */
  plhs[0] = mxCreateStructMatrix(1, 0, 6, nodefields);
  yxml_init(x, XML_STACK, sizeof(XML_STACK));
  
  /* Load XML document */
  if ((nrhs != 1) || (!mxIsChar(prhs[0])))
    mexErrMsgTxt("First input must be a string");
  xml = get_data(prhs[0]);
  
  /* Parse XML document */
  for(doc = xml; *doc; doc++) {
    
    yxml_ret_t r = yxml_parse(x, *doc);
    switch (r) {
      case YXML_OK:
        continue;
      case YXML_CONTENT:
        if (XML_START == NULL && doc[0] != 0x20 && doc[0] != 0x09 && doc[0] != 0x0a && doc[0] != 0x0d)
          XML_START = XML_CURR = doc;
        else
          XML_CURR = doc;
        /* see x->data */
        break;
      case YXML_ELEMSTART:
      case YXML_ELEMEND:
        if (XML_START) {
          resize_struct(plhs[0],++XML_UID);
          mxSetFieldByNumber(plhs[0],XML_UID-1,0,mxCreateString("chardata"));
          XML_CURR[1] = '\0';
          mxSetFieldByNumber(plhs[0],XML_UID-1,1,mxCreateString(XML_START));
          XML_START = NULL;
          mxSetFieldByNumber(plhs[0],XML_UID-1,4,mxCreateDoubleScalar(XML_UID));
          if (XML_CURR_NODE) {
            mxSetFieldByNumber(plhs[0],XML_UID-1,5,mxCreateDoubleScalar(XML_CURR_NODE));
            concatenate_array(plhs[0], XML_CURR_NODE-1, 3, XML_UID);
          }
        }
        switch (r) {
          case YXML_ELEMSTART:
            resize_struct(plhs[0],++XML_UID);
            mxSetFieldByNumber(plhs[0],XML_UID-1,0,mxCreateString("element"));
            mxSetFieldByNumber(plhs[0],XML_UID-1,1,mxCreateString(x->elem));
            mxSetFieldByNumber(plhs[0],XML_UID-1,2,mxCreateStructMatrix(1,0,2,attrfields));
            mxSetFieldByNumber(plhs[0],XML_UID-1,4,mxCreateDoubleScalar(XML_UID));
            if (XML_CURR_NODE) {
              mxSetFieldByNumber(plhs[0],XML_UID-1,5,mxCreateDoubleScalar(XML_CURR_NODE));
              concatenate_array(plhs[0], XML_CURR_NODE-1, 3, XML_UID);
            }
            XML_CURR_NODE = XML_UID;
            break;
          case YXML_ELEMEND:
            mx = mxGetFieldByNumber(plhs[0],XML_CURR_NODE-1,5);
            XML_CURR_NODE = (mx)? (int)mxGetScalar(mx) : 0;
            break;
          default:
            continue;
        }
        break;
      case YXML_ATTRSTART:
        nattr = mxGetNumberOfElements(mxGetFieldByNumber(plhs[0],XML_UID-1,2));
        resize_struct(mxGetFieldByNumber(plhs[0],XML_UID-1,2),++nattr);
        mxSetFieldByNumber(mxGetFieldByNumber(plhs[0],XML_UID-1,2),nattr-1,0,mxCreateString(x->attr));
        break;
      case YXML_ATTRVAL:
        if (XML_START == NULL) XML_START = doc;
        /* see x->data */
        break;
      case YXML_ATTREND:
        nattr = mxGetNumberOfElements(mxGetFieldByNumber(plhs[0],XML_UID-1,2));
        doc[0] = '\0';
        mxSetFieldByNumber(mxGetFieldByNumber(plhs[0],XML_UID-1,2),nattr-1,1,mxCreateString(XML_START));
        XML_START = NULL;
        break;
      case YXML_PISTART:
      case YXML_PICONTENT:
      case YXML_PIEND:
        mexPrintf("Processing Instructions are ignored\n");
        break;
      case YXML_EEOF:
        mexErrMsgTxt("Unexpected EOF");
        break;
      case YXML_EREF:
        mexErrMsgTxt("Invalid character or entity reference");
        break;
      case YXML_ECLOSE:
        mexErrMsgTxt("Close tag does not match open tag");
        break;
      case YXML_ESTACK:
        mexErrMsgTxt("Stack overflow");
        break;
      case YXML_ESYN:
        mexErrMsgTxt("Syntax error");
        break;
      default:
        mexErrMsgTxt("Unexpected error or XML element");
    }
  }
  if (yxml_eof(x) != YXML_OK)
    mexErrMsgTxt("The XML document was not valid");
  
  mxFree(xml);
}
