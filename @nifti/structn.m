function t = structn(obj)
% Convert a NIFTI-1 object into a form of struct
% _______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

%
% $Id: structn.m 174 2005-05-24 11:03:32Z john $


if numel(obj)~=1,
    error('Too many elements to convert');
end;
fn = fieldnames(obj);
for i=1:length(fn)
    tmp = subsref(obj,struct('type','.','subs',fn{i}));
    if ~isempty(tmp)
        t.(fn{i}) = tmp;
    end;
end;
return;
