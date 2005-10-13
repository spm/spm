function disp(obj)
% Disp a NIFTI-1 object
% _______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

%
% $Id: disp.m 253 2005-10-13 15:31:34Z guillaume $


sz = size(obj);
fprintf('NIFTI object: ');
if length(sz)>4,
    fprintf('%d-D\n',length(sz));
else
    for i=1:(length(sz)-1),
        fprintf('%d-by-',sz(i));
    end;
    fprintf('%d\n',sz(end));
end;
if prod(sz)==1,
    display(structn(obj))
end;
return;
