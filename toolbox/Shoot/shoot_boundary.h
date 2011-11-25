
#ifdef NEUMANN
    /* Neumann boundary condition */
    static mwSignedIndex neumann(mwSignedIndex i, mwSize m)
    {
        if (m==1)
            return(0);
        else
        {
            mwSignedIndex m2 = m*2;
            i = (i<0) ? m2-((-i-1)%m2)-1 : (i%m2);
            if (m<=i)
                return(m2-i-1);
            else
                return(i);
        }
    }
#   define BOUND(i,m) neumann(i,m)
#else
    /* circulant boundary condition */
#   define BOUND(i,m) (((signed)(i)>=0) ? (signed)(i)%((signed)m) : (((signed)m)+(signed)(i)%((signed)m))%(signed)m)
#endif

