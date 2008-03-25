function res = reject(obj, ind, flag)
% Method for getting/setting rejection flags
% FORMAT res = reject(obj, ind, flag)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id$

switch nargin
    case 1
        res = getreject(this);        
    case 2
        res = getreject(this, ind)
    case 3
        res = setreject(this, ind, flag);
    otherwise
end

function res = getreject(this, ind)
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

function this = setreject(this, ind, flag)

% include check for datasource!
if obj.Nsamples>0
    if ~isfield(this.trials(1), 'reject')
        % initialize
        [this.trials.reject] = deal(0);
    end
    
    [this.trials(ind).reject] = deal(flag);
end