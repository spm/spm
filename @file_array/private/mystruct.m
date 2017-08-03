function t = mystruct(obj)
%__________________________________________________________________________
% Copyright (C) 2005-2017 Wellcome Trust Centre for Neuroimaging

%
% $Id: mystruct.m 7147 2017-08-03 14:07:01Z spm $


if numel(obj)~=1
    error('Too many elements to convert.');
end
fn = fieldnames(obj);
for i=1:length(fn)
    t.(fn{i}) = subsref(obj,struct('type','.','subs',fn{i}));
end
