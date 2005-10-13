function en = end(a,k,n)
% Overloaded end function for spm_file_array objects.
% _______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id: end.m 253 2005-10-13 15:31:34Z guillaume $

if n>length(a.dim),
	a.dim = [a.dim(1:(k-1)) prod(a.dim(k:end))];
end;
en = a.dim(k);
