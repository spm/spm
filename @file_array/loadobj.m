function b = loadobj(a)
% loadobj for file_array class
%__________________________________________________________________________
% Copyright (C) 2005-2017 Wellcome Trust Centre for Neuroimaging

%
% $Id: loadobj.m 7147 2017-08-03 14:07:01Z spm $


if isa(a,'file_array')
    b = a;
else
    a = permission(a, 'rw');
    b = file_array(a);
end
