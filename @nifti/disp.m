function disp(obj)
% Disp a NIFTI-1 object
%__________________________________________________________________________
% Copyright (C) 2005-2017 Wellcome Trust Centre for Neuroimaging

%
% $Id: disp.m 7147 2017-08-03 14:07:01Z spm $


sz = size(obj);
fprintf('NIFTI object: ');
if length(sz)>4
    fprintf('%d-D\n',length(sz));
else
    for i=1:(length(sz)-1)
        fprintf('%d-by-',sz(i));
    end
    fprintf('%d\n',sz(end));
end
if prod(sz)==1
    disp(structn(obj))
end
