function en = end(a,k,n)
% Overloaded end function for file_array objects.
% _______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

%
% $Id: end.m 253 2005-10-13 15:31:34Z guillaume $

dim = size(a);
if k>length(dim)
    en = 1;
else
    if n<length(dim),
	dim = [dim(1:(n-1)) prod(dim(n:end))];
    end;
    en = dim(k);
end;
