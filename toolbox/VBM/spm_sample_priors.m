function [s,ds1,ds2,ds3] = spm_sample_priors(b,x1,x2,x3,bg)
deg = 3;

if nargout<=1,
    s      = spm_bsplins(b,x1,x2,x3,[deg deg deg  0 0 0]);
    msk    = find(~finite(s));
    s(msk) = bg;
else,
    [s,ds1,ds2,ds3] = spm_bsplins(b,x1,x2,x3,[deg deg deg  0 0 0]);
    msk      = find(~finite(s));
    s(msk)   = bg;
    ds1(msk) = 0;
    ds2(msk) = 0;
    ds3(msk) = 0;
end;
