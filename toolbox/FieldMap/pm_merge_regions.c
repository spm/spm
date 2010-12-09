#include "mex.h"
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <float.h>

#ifndef PI
#define PI 3.14159265358979
#endif

/* changed const int * to const void * DMT 8/22/07 */
int ii_cmp(const void * a,
           const void * b);
/* changed const int * to const void * DMT 8/22/07 */
int jj_cmp(const void * a,
           const void * b);

void merge_regions(/* Input */
                   int     *nn,   /* Number of bordering voxels per pair of regions. */
                   double  *pp,   /* Sum of phase differences along border between regions of pair. */
                   int     nc,    /* Length of list of region pairs. */
                   int     *rs,   /* Size (voxels) per region. */
                   int     nr,    /* Number of regions. */
                   /* Output */
                   int     *rw);  /* Multiple of 2pi wrap for each region. */

int *binsearch(int   *list,
               int   *ia,
               int   n,
               int   post,
               int   *length);

int linsearch(int   *list,
              int   *ia,
              int   n,
              int   post);

int my_round(double  iv);

void sillysort(int  *ii_ia,
               int  *jj_ia,
               int  n);

/* 
   Declare arrays of row and column indicies into
   connectogram matrices as module global to allow
   us to sort on them and create index-arrays using
   stdlib function qsort. It's not my finest hour.
*/

static int     *ii=NULL;
static int     *jj=NULL;
static double  *cc=NULL;

/* changed const int * to const void * DMT 8/22/07 */
int ii_cmp(const void * a,
           const void * b)
{
   /*   if (jj[((int) *a)] > jj[((int) *b)]) {return(1);}
     else if (jj[((int) *a)] == jj[((int) *b)]) {return(0);}
     else {return(-1);} Replaced DMT 8/22/07 */
   if (ii[(*(int*)a)] > ii[(*(int*)b)]) {return(1);}
   else if (ii[(*(int*) a)] == ii[(*(int*) b)]) {return(0);}
   else {return(-1);}
}

/* changed const int * to const void * DMT 8/22/07 */
int jj_cmp(const void * a,
           const void * b)
{
     /*   if (jj[((int) *a)] > jj[((int) *b)]) {return(1);}
     else if (jj[((int) *a)] == jj[((int) *b)]) {return(0);}
     else {return(-1);} Replaced DMT 8/22/07 */
   if (jj[(*(int*)a)] > jj[(*(int*)b)]) {return(1);}
   else if (jj[(*(int*) a)] == jj[(*(int*) b)]) {return(0);}
   else {return(-1);}

}

int my_round(double  iv)
{
  if (iv > 0.0) {return((int) (iv+0.5));}
  else if (iv < 0.0) {return((int) (iv-0.5));}
  else return(0);
}

void merge_regions(/* Input */
                   int     *nn,   /* Number of bordering voxels per pair of regions. */
                   double  *pp,   /* Sum of phase differences along border between regions of pair. */
                   int     nc,    /* Length of list of region pairs. */
                   int     *rs,   /* Size (voxels) per region. */
                   int     nr,    /* Number of regions. */
                   /* Output */
                   int     *rw)   /* Multiple of 2pi wrap for each region. */
{
   double     k=0.0;                     /* Average multiple of 2pi phase difference between regions. */  
   double     *cc=NULL;                  /* Array of "mis-wrap-costs" for each connection. */
   double     maxcc = -FLT_MAX;          /* Current maximum cost. */
   int        l=0;                       /* Wrap when merging olbl with nlbl. */
   int        i=0, j=0, mi=0;            /* Some indicies. */
   int        ci=0, ci2=0;               /* Indexes into ii, jj, nn, pp and cc. */
   int        sign1=1, sign2=1;          /* Used to keep track of sign changes when wrapping. */
   int        cnc=0;                     /* Current remaining no of connections. */ 
   int        ocnc=0;                    /* Remaining no. of connections at start of current iteration. */
   int        row=0, col=0;              /* Row and column (two region labels) for current pair to merge. */
   int        olbl=0;                    /* Label numbers for current pair to */
   int        nlbl=0;                    /* merge such that nlbl+olbl->nlbl.  */
   int        olbli_n=0;                 /* Number of instances of olbl in row-list. */
   int        olblj_n=0;                 /* Number of instances of olbl in column-list. */
   int        *olbli=NULL;               /* List of indexes into ii such that ii[olbli[0:olbli_n-1]]==olbl */
   int        *olblj=NULL;               /* List of indexes into jj such that jj[olblj[0:olblj_n-1]]==olbl */
   int        *mlbli_n=NULL;             /* Number of instances of nlbl in row-list. */
   int        *mlblj_n=NULL;             /* Number of instances of nlbl in column-list. */
   int        **mlbli=NULL;              /* List of indexes into ii such that ii[nlbli[0:nlbli_n-1]]==nlbl */
   int        **mlblj=NULL;              /* List of indexes into jj such that jj[nlblj[0:nlblj_n-1]]==nlbl */
   int        *ii_ia=NULL, *jj_ia=NULL;  /* Index arrays for ii and jj. */
   int        *eqlist=NULL;              /* Equivalence list, indicating the region */
                                         /* i+1 is equivalent to eqlist[i]          */

   /*
   Allocate, initialise and sort index-arrays for 
   ii (row-array) and jj (column array) such that
   ii[ii_ia[:]] is sorted in ascending order as
   is jj[jj_ia[:]].
   */
   ii_ia = (int *) mxCalloc(nc,sizeof(int));
   jj_ia = (int *) mxCalloc(nc,sizeof(int));
   for (i=0; i<nc; i++) {ii_ia[i]=i; jj_ia[i]=i;}
   qsort(ii_ia,nc,sizeof(int),ii_cmp);
   qsort(jj_ia,nc,sizeof(int),jj_cmp);

   /*
   Initialise region-wrap and equivalence lists.
   */
   eqlist = (int *) mxCalloc(nr,sizeof(int));
   for (i=0; i<nr; i++) {rw[i] = 0; eqlist[i] = i+1;}

   /*
   Allocate memory for lists of pointers to stretches
   into ii_ia and jj_ia.
   */
   mlbli_n = (int *) mxCalloc(nc,sizeof(int));   
   mlblj_n = (int *) mxCalloc(nc,sizeof(int));   
   mlbli = (int **) mxCalloc(nc,sizeof(mlbli[0]));   
   mlblj = (int **) mxCalloc(nc,sizeof(mlblj[0]));   

   /*
   Make cc array that is used to determine
   what region pair to merge.
   */
   cc = (double *) mxCalloc(nc,sizeof(double));
   for (i=0; i<nc; i++)
   {
      k = -pp[i]/(2.0*PI*nn[i]);
      cc[i] = nn[i]*(0.5-fabs(k-floor(k+0.5)));
   }
 
   /*
   Go through loop as many times as there are regions, each
   time selecting one pair of regions to merge, and updating
   all statistics regarding connections to remaining regions.
   */
   cnc = nc;
   for (i=0; i<nr; i++)
   {
      ocnc = cnc;  /* Current no. of connections at start of loop. */

      /*
      Find index of next pair to merge.
      */
      for (j=0, mi=0, maxcc=0.0; j<nc; j++)
      {
         if (cc[j] > maxcc) {maxcc=cc[j]; mi=j;}
      }

      /*
      Find row and column of next pair to merge.
      */
      row = ii[mi]; 
      col = jj[mi];

      /*
      Determine which label to give the merged region
      and update region size stats.
      */
      nlbl = (rs[row-1] > rs[col-1]) ? row : col;
      olbl = (nlbl == row) ? col : row;
      rs[nlbl-1] = rs[row-1] + rs[col-1];
      rs[olbl-1] = 0;

      /*
      Determine if region should be wrapped before merging.
      */
      l = my_round(-pp[mi]/(2.0*PI*nn[mi]));

      /*
      Set wrapping for region olbl (and equivalent) in list, 
      and update equivalence list.
      */
      if (nlbl == row) {sign1 = -1;}
      else {sign1 = 1;}
      for (j=0; j<nr; j++)
      {
    if (eqlist[j] == olbl) 
        {
           rw[j] += sign1*l;
           eqlist[j] = nlbl;
        }
      }

      /*
      Find pointers to start of stretches of all pairs 
      that contain the old region as one party.
      */
      olbli = binsearch(ii,ii_ia,cnc,olbl,&olbli_n);
      olblj = binsearch(jj,jj_ia,cnc,olbl,&olblj_n);

      /*
      Find pointers to starts of stretches of all pairs
      that contain members that are somewhere paired
      with the "old" region. We have to do that here while
      ii_ia and jj_ia are still sorted.
      This would appear more complicated than searching
      for stretches containing nlbl, but hopefully it will
      mean that we do linear searches (in the loop below)
      in short lists, allowing us to do binary searches
      here (in the long list).
      */
      for (j=0; j<olbli_n; j++)
      {
     ci = olbli[j];
     if (jj[ci] < nlbl) /* Search in row-list. */
     {
        mlbli[j] = binsearch(ii,ii_ia,cnc,jj[ci],&(mlbli_n[j]));
         }
         else if (jj[ci] > nlbl) /* Search in column list. */
     {
        mlbli[j] = binsearch(jj,jj_ia,cnc,jj[ci],&(mlbli_n[j]));
         }
         else /* Paired with nlbl, skip it. */
     {
        mlbli[j] = NULL;
        mlbli_n[j] = 0;
         }
      }
      for (j=0; j<olblj_n; j++)
      {
     ci = olblj[j];
     if (ii[ci] < nlbl) /* Search in row-list. */
     {
        mlblj[j] = binsearch(ii,ii_ia,cnc,ii[ci],&(mlblj_n[j]));
         }
         else if (ii[ci] > nlbl) /* Search in column list. */
     {
        mlblj[j] = binsearch(jj,jj_ia,cnc,ii[ci],&(mlblj_n[j]));
         }
         else /* Paired with nlbl, skip it. */
     {
        mlblj[j] = NULL;
        mlblj_n[j] = 0;
         }
      }

      /*
      For each of these pairs we can consider three possible cases
      1. It is paired with the new label, and should be deleted.
      2. It is paired with a label that is also in a pair with
         the new label. In this case the stats for the pair
         'new label'-label should be updated and 'old label'-label
         be deleted.
      3. It is paired with a label that is NOT in a pair with
         the new label. In this case the stats for 'old label'-label
         should be transferred to a 'new label'-label pair and the 
         'old label'-label should be deleted.
      */
      /*
      First go through all instances where olbl is row index.
      */
      if (olbli)
      {
     sign2 = 1;
         for (j=0; j<olbli_n; j++)
         {
            ci = olbli[j];
        if (jj[ci] == nlbl)
            {
               nn[ci] = 0; /* Delete */
               pp[ci] = 0.0;
               cc[ci] = 0.0;
               ii[ci] = jj[ci] = INT_MAX;
               cnc--;
            }
            /* 
            Check if the label is currently paired up with nlbl.
            Remember that row index (ii) is always smaller than
            column index->we know where to look.
            */
            else if ((jj[ci] < nlbl) && mlbli[j] && ((ci2=linsearch(jj,mlbli[j],mlbli_n[j],nlbl))>-1))
            /*
            We found the label jj[ci] in the row-list matched with nlbl in the column list.
            Since olbl was in the row-list, this means we should sign-reverse sum of
            differences stats.
            */
            {
           if (l) {pp[ci] = pp[ci] + sign1*sign2*l*2.0*PI*nn[ci];}
               nn[ci2] += nn[ci];
               pp[ci2] -= pp[ci]; /* Third sign reversal. */
               k = -pp[ci2]/(2.0*PI*nn[ci2]);
               cc[ci2] = nn[ci2]*(0.5-fabs(k-floor(k+0.5)));
               nn[ci] = 0; /* Delete */
               pp[ci] = 0.0;
               cc[ci] = 0.0;
               ii[ci] = jj[ci] = INT_MAX;
               cnc--;
            }
            else if ((jj[ci] > nlbl) && mlbli[j] && ((ci2=linsearch(ii,mlbli[j],mlbli_n[j],nlbl))>-1))
            /*
            We found label jj[ci] in the column-list matched with nlbl in the row list.
            Since olbl was also in the row-list, this means we should just add to sum of
            differences stats.
            */
            {
           if (l) {pp[ci] = pp[ci] + sign1*sign2*l*2.0*PI*nn[ci];}
               nn[ci2] += nn[ci];
               pp[ci2] += pp[ci];
               k = -pp[ci2]/(2.0*PI*nn[ci2]);
               cc[ci2] = nn[ci2]*(0.5-fabs(k-floor(k+0.5)));
               nn[ci] = 0; /* Delete */
               pp[ci] = 0.0;
               cc[ci] = 0.0;
               ii[ci] = jj[ci] = INT_MAX;
               cnc--;
            }
            else
            /* 
            So, this label has not been mixed up with nlbl before. Well, it is now.
            */
            {
               if (jj[ci] < nlbl)
               {
          if (l) 
                  {
                     pp[ci] = pp[ci] + sign1*sign2*l*2.0*PI*nn[ci];
                     k = -pp[ci]/(2.0*PI*nn[ci]);
                     cc[ci] = nn[ci]*(0.5-fabs(k-floor(k+0.5)));
                  }
                  ii[ci] = jj[ci];
                  jj[ci] = nlbl;
                  pp[ci] *= -1.0;
               }
               else
               {
          if (l) 
                  {
                     pp[ci] = pp[ci] + sign1*sign2*l*2.0*PI*nn[ci];
                     k = -pp[ci]/(2.0*PI*nn[ci]);
                     cc[ci] = nn[ci]*(0.5-fabs(k-floor(k+0.5)));
                  }
                  ii[ci] = nlbl;
               }
            }
         }
      }
      /*
      Now go through list where olbl is column index.
      */
      if (olblj)
      {
     sign2 = -1;
         for (j=0; j<olblj_n; j++)
         {
            ci = olblj[j];
        if (ii[ci] == nlbl)
            {
               nn[ci] = 0; /* Delete */
               pp[ci] = 0.0;
               cc[ci] = 0.0;
               ii[ci] = jj[ci] = INT_MAX;
               cnc--;
            }
            /* 
            Check if the label is currently paired up with nlbl.
            Remember that row index (ii) is always smaller than
            column index->we know where to look.
            */
            else if ((ii[ci] < nlbl) && mlblj[j] && ((ci2=linsearch(jj,mlblj[j],mlblj_n[j],nlbl))>-1))
            /*
            We found the label ii[ci] in the row-list matched with nlbl in the column list.
            Since olbl was in the column-list, this means we should just add sum of
            differences stats.
            */
            {
           if (l) {pp[ci] = pp[ci] + sign1*sign2*l*2.0*PI*nn[ci];}
               nn[ci2] += nn[ci];
               pp[ci2] += pp[ci]; 
               k = -pp[ci2]/(2.0*PI*nn[ci2]);
               cc[ci2] = nn[ci2]*(0.5-fabs(k-floor(k+0.5)));
               nn[ci] = 0; /* Delete */
               pp[ci] = 0.0;
               cc[ci] = 0.0;
               ii[ci] = jj[ci] = INT_MAX;
               cnc--;
            }
            else if ((ii[ci] > nlbl) && mlblj[j] && ((ci2=linsearch(ii,mlblj[j],mlblj_n[j],nlbl))>-1))
            /*
            We found label ii[ci] in the column-list matched with nlbl in the row list.
            Since olbl was also in the column-list, this means we should sign-reverse when
            adding to sum of differences stats.
            */
            {
           if (l) {pp[ci] = pp[ci] + sign1*sign2*l*2.0*PI*nn[ci];}
               nn[ci2] += nn[ci];
               pp[ci2] -= pp[ci]; /* Third sign reversal. */
               k = -pp[ci2]/(2.0*PI*nn[ci2]);
               cc[ci2] = nn[ci2]*(0.5-fabs(k-floor(k+0.5)));
               nn[ci] = 0; /* Delete */
               pp[ci] = 0.0;
               cc[ci] = 0.0;
               ii[ci] = jj[ci] = INT_MAX;
               cnc--;
            }
            else
            /* 
            So, this label has not been mixed up with nlbl before. Well, it is now.
            */
            {
               if (ii[ci] < nlbl)
               {
          if (l) 
                  {
                     pp[ci] = pp[ci] + sign1*sign2*l*2.0*PI*nn[ci];
                     k = -pp[ci]/(2.0*PI*nn[ci]);
                     cc[ci] = nn[ci]*(0.5-fabs(k-floor(k+0.5)));
                  }
                  jj[ci] = nlbl;
               }
               else
               {
          if (l) 
                  {
                     pp[ci] = pp[ci] + sign1*sign2*l*2.0*PI*nn[ci];
                     k = -pp[ci]/(2.0*PI*nn[ci]);
                     cc[ci] = nn[ci]*(0.5-fabs(k-floor(k+0.5)));
                  }
                  jj[ci] = ii[ci];
                  ii[ci] = nlbl;
                  pp[ci] *= -1.0;
               }
            }
         }
      }
      if (!cnc) {break;} /* No more connections, we're done. */
      /*
      Resort index arrays ii_ia and jj_ia. Connections that are
      no longer valid have been given ii and jj values of INT_MAX,
      which means they will be out at the end of the ii_ia and jj_ia
      arrays, and will not be considered in subsequent iterations.
      */
      if ((ocnc-cnc) == 1)
      /*
      A large proportion of mergings involve merging of regions
      that have connections only with the current nlbl. Clearly
      we can do much better than a full quicksort for those cases.
      */
      {
     sillysort(ii_ia,jj_ia,ocnc);
      }
      else
      {
         qsort(ii_ia,ocnc,sizeof(int),ii_cmp);
         qsort(jj_ia,ocnc,sizeof(int),jj_cmp);
      }

      /*      return; */
   }

   mxFree(cc);
   mxFree(ii_ia);
   mxFree(jj_ia);   
   mxFree(mlbli_n);
   mxFree(mlblj_n);
   mxFree(mlbli);
   mxFree(mlblj);
   return;
}


/*
** Quick and dirty routine to sort lists when there is only one
** (meaning that there is one INT_MAX wrong sorted in each
** index array) element sorted wrong. 
**It is all very primitive. A slightly more informative name 
** would be "supersillysort", though a tad bit long.
**/

void sillysort(int  *ii_ia,
               int  *jj_ia,
               int  n)
{
   int    i=0, j=0;
   int    tmp=0;

   for (i=0; i<n; i++)
   {
      if (ii[ii_ia[i]] == INT_MAX)
      {
     tmp = ii_ia[i];
         for (j=i; j<(n-1); j++)
     {
        ii_ia[j] = ii_ia[j+1];
         }
         ii_ia[n-1] = tmp;
         break;
      }
   }

   for (i=0; i<n; i++)
   {
      if (jj[jj_ia[i]] == INT_MAX)
      {
     tmp = jj_ia[i];
         for (j=i; j<(n-1); j++)
     {
        jj_ia[j] = jj_ia[j+1];
         }
         jj_ia[n-1] = tmp;
         break;
      }
   }

   return;
}

/*
** Routine to do binary search in an integer list sorted
** in ascending order by an index-array. Returns pointer into
** index-array for the lowest i for which list[ia[i]]==post. 
** A little more convenient than using bsearch
** considering the complication of the index-array.
** It is used to find base and length of stretches where
** list[ia[i]]==post. 
*/ 
int *binsearch(int   *list,
               int   *ia,
               int   n,
               int   post,
               int   *length)
{
   int   ll=-1, ul=n;
   int   mp=0;

   while (ul-ll > 1)
   {
      mp = ll + (ul-ll) / 2;
      if (list[ia[mp]] == post) {ll=mp; break;}
      else if (list[ia[mp]] < post) {ll=mp;}
      else {ul=mp;}
   }
   if ((ll<0) || (list[ia[ll]] != post)) {(*length)=0; return(NULL);}
   else
   {
      while (ll && (list[ia[ll-1]] == post)) {ll--;}
   }
   (*length)=0; 
   while (((ll+(*length)) < n) && list[ia[ll+(*length)]] == post) {(*length)++;}

   return(&(ia[ll]));
}

/*
**
** Routine to do linear search for single instance of item "post".
**
*/ 

int linsearch(int   *list,
              int   *ia,
              int   n,
              int   post)
{
   int   i=0;

   for (i=0; i<n; i++)
   {
      if (list[ia[i]] == post) {return(ia[i]); }
   }

   return(-1);
}



/* Gateway function with error check. */

void mexFunction(int             nlhs,      /* No. of output arguments */
                 mxArray         *plhs[],   /* Output arguments. */ 
                 int             nrhs,      /* No. of input arguments. */
                 const mxArray   *prhs[])   /* Input arguments. */
{
   int            i=0;
   int            ndim=0, lima_ndim=0;
   int            nc=0, nr=0, ni=0;
   const int      *cdim=NULL, *lima_cdim=NULL;
   unsigned int   dim[3];
   int            *nn=NULL;
   int            *rs=NULL;
   int            *rw=NULL;
   int            *lima=NULL;
   double         *pp=NULL;
   double         *pm=NULL;
   double         *d_tmp=NULL;
   double         *opm=NULL;

   if (nrhs == 0) mexErrMsgTxt("usage: pm = pm_merge_regions(pm,lima,i,j,n,p,rs)");
   if (nrhs != 7) mexErrMsgTxt("pm_merge_regions: 7 input arguments required");
   /*
   if (nlhs != 1) mexErrMsgTxt("pm_merge_regions: 1 output argument required");
   */

   /* Get phase map. */

   if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]))
   {
      mexErrMsgTxt("pm_merge_regions: pm must be numeric, real, full and double");
   }
   ndim = mxGetNumberOfDimensions(prhs[0]);
   if ((ndim < 2) | (ndim > 3))
   {
      mexErrMsgTxt("pm_merge_regions: pm must be 2 or 3-dimensional");
   }
   cdim = mxGetDimensions(prhs[0]);
   pm = mxGetPr(prhs[0]);

   /* Get image of labels. */

   if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1]))
   {
      mexErrMsgTxt("pm_merge_regions: lima must be numeric, real, full and double");
   }
   lima_ndim = mxGetNumberOfDimensions(prhs[1]);
   if (lima_ndim != ndim)
   {
      mexErrMsgTxt("pm_merge_regions: pm and lima must have same dimensionality");
   }
   lima_cdim = mxGetDimensions(prhs[1]);
   for (i=0; i<ndim; i++)
   {
      if (cdim[i] != lima_cdim[i])
      {
         mexErrMsgTxt("pm_merge_regions: pm and lima must have same size");
      }
   }
   d_tmp = mxGetPr(prhs[1]);

   /* Fix dimensions to allow for 2D and 3D data. */

   dim[0]=cdim[0]; dim[1]=cdim[1];
   if (ndim==2) {dim[2]=1; ndim=3;} else {dim[2]=cdim[2];} 
   for (i=0, ni=1; i<ndim; i++)
   {
      ni *= dim[i];
   }

   /* Convert double representation of lima into int's. */

   lima = (int *) mxCalloc(ni,sizeof(int));
   for (i=0; i<ni; i++) {lima[i] = ((int) (d_tmp[i]+0.1));} 

   /* Get vector of row indicies into connectogram matrix */

   if (!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2]) || mxIsSparse(prhs[2]) || !mxIsDouble(prhs[2]))
   {
      mexErrMsgTxt("pm_merge_regions: i must be numeric, real, full and double");
   }
   nc = mxGetM(prhs[2]);
   d_tmp = mxGetPr(prhs[2]);
   if (mxGetN(prhs[2]) != 1)
   {
      mexErrMsgTxt("pm_merge_regions: i must be a column matrix");
   }
   ii = (int *) mxCalloc(nc,sizeof(int));
   for (i=0; i<nc; i++) {ii[i] = ((int) (d_tmp[i]+0.1));}

   /* Get vector of column indicies into connectogram matrix */

   if (!mxIsNumeric(prhs[3]) || mxIsComplex(prhs[3]) || mxIsSparse(prhs[3]) || !mxIsDouble(prhs[3]))
   {
      mexErrMsgTxt("pm_merge_regions: j must be numeric, real, full and double");
   }
   if (mxGetM(prhs[3]) != nc || mxGetN(prhs[3]) != 1)
   {
      mexErrMsgTxt("pm_merge_regions: j must be a column matrix of same size as i");
   }
   d_tmp = mxGetPr(prhs[3]);
   jj = (int *) mxCalloc(nc,sizeof(int));
   for (i=0; i<nc; i++) {jj[i] = ((int) (d_tmp[i]+0.1));}

   /* Get entrys for first connectogram matrix (containing border sizes). */

   if (!mxIsNumeric(prhs[4]) || mxIsComplex(prhs[4]) || mxIsSparse(prhs[4]) || !mxIsDouble(prhs[4]))
   {
      mexErrMsgTxt("pm_merge_regions: n must be numeric, real, full and double");
   }
   if (mxGetM(prhs[4]) != nc || mxGetN(prhs[4]) != 1)
   {
      mexErrMsgTxt("pm_merge_regions: n must be a column matrix of same size as i");
   }
   d_tmp = mxGetPr(prhs[4]);
   nn = (int *) mxCalloc(nc,sizeof(int));
   for (i=0; i<nc; i++) {nn[i] = ((int) (d_tmp[i]+0.1));}

   /* Get entrys for second connectogram matrix (containing phase differences along borders). */

   if (!mxIsNumeric(prhs[5]) || mxIsComplex(prhs[5]) || mxIsSparse(prhs[5]) || !mxIsDouble(prhs[5]))
   {
      mexErrMsgTxt("pm_merge_regions: p must be numeric, real, full and double");
   }
   if (mxGetM(prhs[5]) != nc || mxGetN(prhs[5]) != 1)
   {
      mexErrMsgTxt("pm_merge_regions: p must be a column matrix of same size as i");
   }
   d_tmp = mxGetPr(prhs[5]);
   pp = (double *) mxCalloc(nc,sizeof(double)); /* Use local copy to avoid side effects. */
   memcpy(pp,d_tmp,nc*sizeof(double));   

   /* Get vector of region sizes (in voxels) */

   if (!mxIsNumeric(prhs[6]) || mxIsComplex(prhs[6]) || mxIsSparse(prhs[6]) || !mxIsDouble(prhs[6]))
   {
      mexErrMsgTxt("pm_merge_regions: rs must be numeric, real, full and double");
   }
   nr = mxGetM(prhs[6]);
   d_tmp = mxGetPr(prhs[6]);
   if (mxGetN(prhs[6]) != 1)
   {
      mexErrMsgTxt("pm_merge_regions: rs must be a column matrix");
   }
   rs = (int *) mxCalloc(nr,sizeof(int));
   for (i=0; i<nr; i++) {rs[i] = ((int) (d_tmp[i]+0.1));}
   
   /* Allocate mem for unwrapped output phasemap. */

      plhs[0] = mxCreateNumericArray(mxGetNumberOfDimensions(prhs[0]),
                                  mxGetDimensions(prhs[0]),mxDOUBLE_CLASS,mxREAL);
   opm = mxGetPr(plhs[0]);
   
   rw = (int *) mxCalloc(nr,sizeof(int));

   merge_regions(nn,pp,nc,rs,nr,rw);

   for (i=0; i<ni; i++)
   {
      if (lima[i]) {opm[i] = pm[i] + 2.0*PI*rw[lima[i]-1];}
      else {opm[i] = pm[i];}
   }

   mxFree(lima);     
   mxFree(ii);     
   mxFree(jj);     
   mxFree(nn);
   mxFree(pp);     
   mxFree(rs);     

   return;
}



