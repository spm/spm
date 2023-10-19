/*
 * John Ashburner
 * Copyright (C) 2023 Wellcome Centre for Human Neuroimaging
 */

/* Cholesky decomposition
 * n  - dimension of matrix a
 * a  - an n \times n matrix
 * p  - an n \times 1 vector
 *
 * A triangle of the input matrix is partially overwritten
 * by the output. Diagonal elements are stored in p.
 */
__device__ void choldcf(USIZE_t n, float a[], /*@out@*/ float p[])
{
    USIZE_t i, j;
    SSIZE_t k;

    for(i=0; i<n; i++)
    {
        float *ai = a + i*n;
        for(j=i; j<n; j++)
        {
            float  t  = ai[j];
            float *aj = a + j*n;
            for(k=(SSIZE_t)i-1; k>=0; k--)
               t -= ai[k]*aj[k];
            if(j==i)
                p[i] = sqrtf(t);
            else
                aj[i] = t/p[i];
        }
    }
}


/* Solve a least squares problem with the results from a
 * Cholesky decomposition
 *
 * n     - Dimension of matrix and data.
 * a & p - Cholesky decomposed matrix.
 * b     - Vector of input data.
 * x     - Vector or outputs.
 */
__device__ void chollsf(USIZE_t n, const float a[], const float p[],
                        const float b[], /*@out@*/ float x[])
{
    SSIZE_t i, k;

    for(i=0; i<(SSIZE_t)n; i++)
    {
        float t = b[i];
        for(k=i-1; k>=0; k--)
            t -= a[i*n+k]*x[k];
        x[i] = t/p[i];
    }
    for(i=(SSIZE_t)n-1; i>=0; i--)
    {
        float t = x[i];
        for(k=i+1; k<(SSIZE_t)n; k++)
            t -= a[k*n+i]*x[k];
        x[i] = t/p[i];
    }
}
