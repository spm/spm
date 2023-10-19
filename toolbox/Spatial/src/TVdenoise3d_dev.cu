/*
 * John Ashburner
 * Copyright (C) 2023 Wellcome Centre for Human Neuroimaging
 */

#define SQUARE(x) (_t=(x), _t*_t)
#define SMALL 1e-8f 

__device__ void TVdenoise3d_fast_dev(USIZE_t i, USIZE_t j, USIZE_t k, float y[], const float x[], const USIZE_t d[],
                                     const float vox[3], const float lambdap[], const float lambdal[])
{
    SSIZE_t d0 = (SSIZE_t)d[0], d1 = (SSIZE_t)d[1], d01 = d0*d1;
    USIZE_t ijk = i + d0*(j + d1*k), n = d0*d1*d[2], m;
    float   _t, ws, w0 = 0.0f, w1 = 0.0f, w2 = 0.0f, w3 = 0.0f;
    float   *yp, yb[MAXVOL][6], *ybp;
    float rv0 = 1.0f/vox[0]/vox[0], rv1 = 1.0f/vox[1]/vox[1], rv2 = 1.0f/vox[2]/vox[2];

    for(m=0, yp=y+ijk; m<d[3]; m++, yp+=n)
    {
        float lambdap_m = lambdap[m];
        float yt[3][3][3];

        yt[0][0][0] = yp[-1-d0-d01]; yt[0][0][1] = yp[  -d0-d01]; yt[0][0][2] = yp[ 1-d0-d01];
        yt[0][1][0] = yp[-1   -d01]; yt[0][1][1] = yp[     -d01]; yt[0][1][2] = yp[ 1   -d01];
        yt[0][2][0] = yp[-1+d0-d01]; yt[0][2][1] = yp[  +d0-d01]; yt[0][2][2] = yp[ 1+d0-d01];

        yt[1][0][0] = yp[-1-d0    ]; yt[1][0][1] = yp[  -d0    ]; yt[1][0][2] = yp[ 1-d0    ];
        yt[1][1][0] = yp[-1       ]; yt[1][1][1] = *yp          ; yt[1][1][2] = yp[ 1       ];
        yt[1][2][0] = yp[-1+d0    ]; yt[1][2][1] = yp[   d0    ]; yt[1][2][2] = yp[ 1+d0    ];

        yt[2][0][0] = yp[-1-d0+d01]; yt[2][0][1] = yp[  -d0+d01]; yt[2][0][2] = yp[ 1-d0+d01];
        yt[2][1][0] = yp[-1   +d01]; yt[2][1][1] = yp[      d01]; yt[2][1][2] = yp[ 1   +d01];
        yt[2][2][0] = yp[-1+d0+d01]; yt[2][2][1] = yp[  +d0+d01]; yt[2][2][2] = yp[ 1+d0+d01];

        ybp    = yb[m];
        ybp[0] = yt[0][1][1];
        ybp[1] = yt[1][0][1];
        ybp[2] = yt[1][1][0];
        ybp[3] = yt[1][1][2];
        ybp[4] = yt[1][2][1];
        ybp[5] = yt[2][1][1];

        w0   += lambdap_m*(SQUARE(yt[1][1][1]-yt[0][1][1])*rv2 + SQUARE(yt[1][1][1]-yt[1][0][1])*rv1 + SQUARE(yt[1][1][1]-yt[1][1][0])*rv0 + SMALL);
        w1   += lambdap_m*(SQUARE(yt[1][1][1]-yt[1][1][2])*rv0 + SQUARE(yt[1][1][2]-yt[1][0][2])*rv1 + SQUARE(yt[1][1][2]-yt[0][1][2])*rv2 + SMALL);
        w2   += lambdap_m*(SQUARE(yt[1][1][1]-yt[1][2][1])*rv1 + SQUARE(yt[1][2][1]-yt[1][2][0])*rv0 + SQUARE(yt[1][2][1]-yt[0][2][1])*rv2 + SMALL);
        w3   += lambdap_m*(SQUARE(yt[1][1][1]-yt[2][1][1])*rv2 + SQUARE(yt[2][1][1]-yt[2][0][1])*rv1 + SQUARE(yt[2][1][1]-yt[2][1][0])*rv0 + SMALL);
    }

    /* See https://francisbach.com/the-%ce%b7-trick-or-the-effectiveness-of-reweighted-least-squares */
    w0 = 1.0f/sqrt(w0);
    w1 =  rv0/sqrt(w1);
    w2 =  rv1/sqrt(w2);
    w3 =  rv2/sqrt(w3);
    ws = w0*(rv2 + rv1 + rv0) + w1 + w2 + w3;

    x += ijk;
    for(m=0, yp=y+ijk; m<d[3]; m++, x+=n, yp+=n)
    {
        ybp = yb[m];
        *yp = (lambdap[m]*((ybp[0]*rv2 + ybp[1]*rv1 + ybp[2]*rv0)*w0 + ybp[3]*w1 + ybp[4]*w2 + ybp[5]*w3) + *x*lambdal[m])
             /(lambdap[m]*ws + lambdal[m]);
    }
}




__device__ void TVdenoise3d_dev(USIZE_t i, USIZE_t j, USIZE_t k, float y[], const float x[], const USIZE_t d[],
                                const float vox[3], const float lambdap[], const float lambdal[])
{
    SSIZE_t d0 = (SSIZE_t)d[0], d1 = (SSIZE_t)d[1], d01 = d0*d1;
    USIZE_t ijk = i + d0*(j + d1*k), n = d0*d1*d[2], m;
    float   w111,  w011, w211, w101, w121, w110, w112;
    float   *yp, yb[MAXVOL][6], *ybp;
    float eta[32];
    float rv0 = 1.0f/vox[0]/vox[0], rv1 = 1.0f/vox[1]/vox[1], rv2 = 1.0f/vox[2]/vox[2];
    /*
        Solving the L1 regularisation problem involves minimising |y| = \frac{y^2}{2\eta} + \frac{\eta}{2}.
        This is achieved by alternating between:
            \hat{\eta} = |\hat{y}|
        and:
            \hat{y} = \argmin_y \frac{y^2}{2\hat{\eta}}

       Prblem here involves multi-channel TV, so instead of |y|, we could minimise
           \sqrt{\sum_m \lambdap_m ((y_{i,j,k,m}-y_{i+1,j,k,m})^2 + (y_{i,j,k,m}-y_{i,j+1,k,m})^2 + (y_{i,j,k,m}-y_{i,j,k+1,m})^2))}

       Also note that the above configuration could involve eight different permutations of neighbours.
       If the above neighbourhood configuration is denoted by +++, we also have ++-, +-+, +--, -++, -+-, ---, --+.
       Therefore, in practice this code uses the average of the eight possible neighbourhood configurations.
     */
    for(m=0; m<32; m++) eta[m] = 0.0f;

    for(m=0, yp=y+ijk; m<d[3]; m++, yp+=n)
    {
        float lambdap_m = lambdap[m];
        float yt[3][3][3];
        float yc, _t, t, s, t0a, t0b, t1a, t1b, t2a, t2b;

        /* Get the 3x3x3 patch */
        yt[0][0][0] = yp[-1-d0-d01]; yt[0][0][1] = yp[  -d0-d01]; yt[0][0][2] = yp[ 1-d0-d01];
        yt[0][1][0] = yp[-1   -d01]; yt[0][1][1] = yp[     -d01]; yt[0][1][2] = yp[ 1   -d01];
        yt[0][2][0] = yp[-1+d0-d01]; yt[0][2][1] = yp[  +d0-d01]; yt[0][2][2] = yp[ 1+d0-d01];

        yt[1][0][0] = yp[-1-d0    ]; yt[1][0][1] = yp[  -d0    ]; yt[1][0][2] = yp[ 1-d0    ];
        yt[1][1][0] = yp[-1       ]; yt[1][1][1] = *yp          ; yt[1][1][2] = yp[ 1       ];
        yt[1][2][0] = yp[-1+d0    ]; yt[1][2][1] = yp[   d0    ]; yt[1][2][2] = yp[ 1+d0    ];

        yt[2][0][0] = yp[-1-d0+d01]; yt[2][0][1] = yp[  -d0+d01]; yt[2][0][2] = yp[ 1-d0+d01];
        yt[2][1][0] = yp[-1   +d01]; yt[2][1][1] = yp[      d01]; yt[2][1][2] = yp[ 1   +d01];
        yt[2][2][0] = yp[-1+d0+d01]; yt[2][2][1] = yp[  +d0+d01]; yt[2][2][2] = yp[ 1+d0+d01];

        yc     = *yp;         /* Central voxel */
        ybp    = yb[m];       /* Pointer for this volume */
        ybp[0] = yt[0][1][1]; /* -z */
        ybp[1] = yt[2][1][1]; /* +z */
        ybp[2] = yt[1][0][1]; /* -y */
        ybp[3] = yt[1][2][1]; /* +y */
        ybp[4] = yt[1][1][0]; /* -x */
        ybp[5] = yt[1][1][2]; /* +x */

        /* Immediate six neighbours */
        t2a    = SQUARE(yc-ybp[0])*rv2;
        t2b    = SQUARE(yc-ybp[1])*rv2;
        t1a    = SQUARE(yc-ybp[2])*rv1;
        t1b    = SQUARE(yc-ybp[3])*rv1;
        t0a    = SQUARE(yc-ybp[4])*rv0;
        t0b    = SQUARE(yc-ybp[5])*rv0;

        /* Centre voxel connecting to three neighbours (eight permutations). */
        s      = t1a + t2a + SMALL; eta[0] += lambdap_m*(s + t0a); eta[1] += lambdap_m*(s + t0b);
        s      = t1b + t2a + SMALL; eta[2] += lambdap_m*(s + t0a); eta[3] += lambdap_m*(s + t0b);
        s      = t1a + t2b + SMALL; eta[4] += lambdap_m*(s + t0a); eta[5] += lambdap_m*(s + t0b);
        s      = t1b + t2b + SMALL; eta[6] += lambdap_m*(s + t0a); eta[7] += lambdap_m*(s + t0b);

        /* Centre voxel connecting to one neighbour, but that neighbour connects to two other voxels
           (four permutations for each connecting neighbour). */
        t        = ybp[0];  /* -z neighbour */
        s        = SMALL  + t2a + SQUARE(t-yt[0][0][1])*rv1;
        eta[8]  += lambdap_m*(s + SQUARE(t-yt[0][1][0])*rv0);
        eta[9]  += lambdap_m*(s + SQUARE(t-yt[0][1][2])*rv0);
        s        = SMALL  + t2a + SQUARE(t-yt[0][2][1])*rv1;
        eta[10] += lambdap_m*(s + SQUARE(t-yt[0][1][0])*rv0);
        eta[11] += lambdap_m*(s + SQUARE(t-yt[0][1][2])*rv0);

        t        = ybp[1];  /* +z */
        s        = SMALL  + t2b + SQUARE(t-yt[2][1][0])*rv1;
        eta[12] += lambdap_m*(s + SQUARE(t-yt[2][0][1])*rv0);
        eta[13] += lambdap_m*(s + SQUARE(t-yt[2][2][1])*rv0);
        s        = SMALL  + t2b + SQUARE(t-yt[2][1][2])*rv1;
        eta[14] += lambdap_m*(s + SQUARE(t-yt[2][0][1])*rv0);
        eta[15] += lambdap_m*(s + SQUARE(t-yt[2][2][1])*rv0);

        t        = ybp[2];  /* -y */
        s        = SMALL  + t1a + SQUARE(t-yt[0][0][1])*rv2;
        eta[16] += lambdap_m*(s + SQUARE(t-yt[1][0][0])*rv0);
        eta[17] += lambdap_m*(s + SQUARE(t-yt[1][0][2])*rv0);
        s        = SMALL  + t1a + SQUARE(t-yt[2][0][1])*rv2;
        eta[18] += lambdap_m*(s + SQUARE(t-yt[1][0][0])*rv0);
        eta[19] += lambdap_m*(s + SQUARE(t-yt[1][0][2])*rv0);

        t        = ybp[3];  /* +y */
        s        = SMALL  + t1b + SQUARE(t-yt[0][2][1])*rv2;
        eta[20] += lambdap_m*(s + SQUARE(t-yt[1][2][0])*rv0);
        eta[21] += lambdap_m*(s + SQUARE(t-yt[1][2][2])*rv0);
        s        = SMALL  + t1b + SQUARE(t-yt[2][2][1])*rv2;
        eta[22] += lambdap_m*(s + SQUARE(t-yt[1][2][0])*rv0);
        eta[23] += lambdap_m*(s + SQUARE(t-yt[1][2][2])*rv0);

        t        = ybp[4];  /* -x */
        s        = SMALL +  t0a + SQUARE(t-yt[0][1][0])*rv2;
        eta[24] += lambdap_m*(s + SQUARE(t-yt[1][0][0])*rv1);
        eta[25] += lambdap_m*(s + SQUARE(t-yt[1][2][0])*rv1);
        s        = SMALL +  t0a + SQUARE(t-yt[2][1][0])*rv2;
        eta[26] += lambdap_m*(s + SQUARE(t-yt[1][0][0])*rv1);
        eta[27] += lambdap_m*(s + SQUARE(t-yt[1][2][0])*rv1);

        t        = ybp[5];  /* +x */
        s        = SMALL +  t0b + SQUARE(t-yt[0][1][2])*rv2;
        eta[28] += lambdap_m*(s + SQUARE(t-yt[1][0][2])*rv1);
        eta[29] += lambdap_m*(s + SQUARE(t-yt[1][2][2])*rv1);
        s        = SMALL +  t0b + SQUARE(t-yt[2][1][2])*rv2;
        eta[30] += lambdap_m*(s + SQUARE(t-yt[1][0][2])*rv1);
        eta[31] += lambdap_m*(s + SQUARE(t-yt[1][2][2])*rv1);
    }

    /* See https://francisbach.com/the-%ce%b7-trick-or-the-effectiveness-of-reweighted-least-squares */
    for(m=0; m<32; m++) eta[m] = 1.0f/sqrt(eta[m]);

    /* Weights from an average of the eight different arrangements of neighbours. */
    w011 = rv2*(eta[ 8]+eta[ 9]+eta[10]+eta[11] + eta[0]+eta[1]+eta[2]+eta[3])*0.125f; /* -z neighbour */
    w211 = rv2*(eta[12]+eta[13]+eta[14]+eta[15] + eta[4]+eta[5]+eta[6]+eta[7])*0.125f; /* +z */
    w101 = rv1*(eta[16]+eta[17]+eta[18]+eta[19] + eta[0]+eta[1]+eta[4]+eta[5])*0.125f; /* -y */
    w121 = rv1*(eta[20]+eta[21]+eta[22]+eta[23] + eta[2]+eta[3]+eta[6]+eta[7])*0.125f; /* +y */
    w110 = rv0*(eta[24]+eta[25]+eta[26]+eta[27] + eta[0]+eta[2]+eta[4]+eta[6])*0.125f; /* -x */
    w112 = rv0*(eta[28]+eta[29]+eta[30]+eta[31] + eta[1]+eta[3]+eta[5]+eta[7])*0.125f; /* +x */
    w111 = w011 + w211 + w101 + w121 + w110 + w112; /* Centre weight */

    x += ijk;
    for(m=0, yp=y+ijk; m<d[3]; m++, x+=n, yp+=n)
    {
        ybp = yb[m];
        *yp = (lambdap[m]*(ybp[0]*w011 + ybp[1]*w211 + ybp[2]*w101 + ybp[3]*w121 + ybp[4]*w110 + ybp[5]*w112) + *x*lambdal[m])
             /(lambdap[m]*w111 + lambdal[m]);
    }
}
