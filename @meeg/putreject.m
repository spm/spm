function obj = putreject(obj, ind, flag)
% Method for chaing a rejection flag
% FORMAT res = putreject(obj, ind, flag)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel

if obj.Nsamples>0
    if ~isfield(obj.trials(1), 'reject')
        [obj.trials.reject] = deal(0);
    end

    [obj.trials(ind).reject] = deal(flag);
end