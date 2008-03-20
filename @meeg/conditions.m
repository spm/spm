function res = conditions(obj, ind)
% Method for getting condition labels, over trials
% FORMAT res = conditions(obj)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
%$Id: conditions.m 1236 2008-03-20 18:15:33Z stefan $

if obj.Nsamples>0
    res = {obj.trials.label};
else
    res = [];
end

if nargin > 1
    res = res(ind);
end

res = strvcat(res);
