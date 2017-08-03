function t = structn(obj)
% Convert a NIFTI-1 object into a form of struct
%__________________________________________________________________________
% Copyright (C) 2005-2017 Wellcome Trust Centre for Neuroimaging

%
% $Id: structn.m 7147 2017-08-03 14:07:01Z spm $


if numel(obj)~=1
    error('Too many elements to convert.');
end
fn = fieldnames(obj);
for i=1:length(fn)
    tmp = subsref(obj,struct('type','.','subs',fn{i}));
    if ~isempty(tmp)
        t.(fn{i}) = tmp;
    end
end
