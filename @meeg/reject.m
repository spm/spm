function res = reject(this, ind, flag)
% Method for getting/setting rejection flags
% FORMAT res = reject(this, ind, flag)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: reject.m 1373 2008-04-11 14:24:03Z spm $

switch nargin
    case 1
        res = getreject(this);        
    case 2
        res = getreject(this, ind);
    case 3
        res = setreject(this, ind, flag);
    otherwise
end

function res = getreject(this, ind)
if this.Nsamples>0
    if isfield(this.trials(1), 'reject')
        res = cat(1,this.trials.reject);
    else
        res = zeros(length(this.trials),1);
    end
else
    res = [];
end

if nargin > 1
    res = res(ind);
end

function this = setreject(this, ind, flag)

% include check for datasource!
if this.Nsamples>0
    if ~isfield(this.trials(1), 'reject')
        % initialize
        [this.trials.reject] = deal(0);
    end
    
    [this.trials(ind).reject] = deal(flag);
end