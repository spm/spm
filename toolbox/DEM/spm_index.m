function [ind] = spm_index(siz,ndx)
% Multiple subscripts from linear index
% FORMAT [ind] = spm_index(siz,ndx)
% siz  - array size
% ndx  - linear index
%
% see: SUB2IND, FIND.
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2012-2022 Wellcome Centre for Human Neuroimaging 

% deal with scalars
%--------------------------------------------------------------------------
k    = cumprod(siz);
if k == 1
    ind = 1;
    return
end

% deal with arrays
%--------------------------------------------------------------------------
len  = length(siz);
if len > 2
    for i   = len:-1:3
        vi  = rem(ndx-1, k(i-1)) + 1;
        vj  = (ndx - vi)/k(i-1) + 1;
        ind(i-2) = vj;
        ndx = vi;
    end
else
    ind = [];
end
if len >= 2
    vi = rem(ndx - 1, siz(1)) + 1;
    v2 = (ndx - vi)/siz(1) + 1;
    v1 = vi;
else 
    v1 = ndx;
end
ind = [v1,v2,ind];


return