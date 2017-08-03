function out = isnan(fa)
% Convert to numeric form
% FORMAT isnan(fa)
% fa - a file_array
%__________________________________________________________________________
% Copyright (C) 2005-2017 Wellcome Trust Centre for Neuroimaging

%
% $Id: isnan.m 7147 2017-08-03 14:07:01Z spm $


bs  = 10240;
m   = size(fa);
fa  = reshape(fa,prod(m),1);
n   = prod(m);
out = false(m);
for i=1:ceil(n/bs)
    ii      = ((((i-1)*bs)+1):min((i*bs),n))';
    tmp     = subsref(fa,struct('type','()','subs',{{ii}}));
    out(ii) = isnan(tmp);
end
