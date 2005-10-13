function t = mystruct(obj)
% _______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

%
% $Id: mystruct.m 253 2005-10-13 15:31:34Z guillaume $


if numel(obj)~=1,
    error('Too many elements to convert');
end;
fn = fieldnames(obj);
for i=1:length(fn)
    t.(fn{i}) = subsref(obj,struct('type','.','subs',fn{i}));
end;
return;
