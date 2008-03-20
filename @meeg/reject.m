function res = reject(obj, ind)
% Method for getting all rejection flags, over trials
% FORMAT res = reject(obj)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel

if obj.Nsamples>0
    if isfield(obj.trials(1), 'reject')
        res = cat(1,obj.trials.reject);
    else
        res = zeros(length(obj.trials),1);
    end
else
    res = [];
end

if nargin > 1
    res = res(ind);
end