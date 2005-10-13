function disp(obj)
% Display a file_array object
% _______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

%
% $Id: disp.m 253 2005-10-13 15:31:34Z guillaume $


if numel(struct(obj))>1,
    fprintf('       %s object: ', class(obj));
    sz = size(obj);
    if length(sz)>4,
        fprintf('%d-D\n',length(sz));
    else
        for i=1:(length(sz)-1),
            fprintf('%d-by-',sz(i));
        end;
        fprintf('%d\n',sz(end));
    end;
else
    display(mystruct(obj))
end;
return;
%=======================================================================

%=======================================================================
function t = mystruct(obj)
fn = fieldnames(obj);
for i=1:length(fn)
    t.(fn{i}) = subsref(obj,struct('type','.','subs',fn{i}));
end;
return;
%=======================================================================
