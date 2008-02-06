function en = end(a,k,n)
% Overloaded end function for file_array objects.
% _______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

%
% $Id: end.m 1137 2008-02-06 15:58:21Z spm $

dim = size(a);
if k>length(dim)
    en = 1;
else
    if n<length(dim),
    dim = [dim(1:(n-1)) prod(dim(n:end))];
    end;
    en = dim(k);
end;
