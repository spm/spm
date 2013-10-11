function [vX] = spm_vec(varargin)
% vectorises a numeric, cell or structure array
% FORMAT [vX] = spm_vec(X);
% X  - numeric, cell or stucture array[s]
% vX - vec(X)
%__________________________________________________________________________
%
% e.g.:
% spm_vec({eye(2) 3}) = [1 0 0 1 3]'
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_vec.m 5691 2013-10-11 16:53:00Z karl $

% initialise X and vX
%--------------------------------------------------------------------------
if nargin == 1
    X = varargin{1};
else
    X = varargin;
end


% vectorise numerical arrays
%--------------------------------------------------------------------------
if isnumeric(X)
    vX = X(:);

% vectorise logical arrays
%--------------------------------------------------------------------------
elseif islogical(X)
    vX = X(:);

% vectorise structure into cell arrays
%--------------------------------------------------------------------------
elseif isstruct(X)
    vX = [];
    f   = fieldnames(X);
    X    = X(:);
    for i = 1:numel(f)
        vX = cat(1,vX,spm_vec({X.(f{i})}));
    end

% vectorise cells into numerical arrays
%--------------------------------------------------------------------------
elseif iscell(X)
    vX   = [];
    for i = 1:numel(X)
        vX = cat(1,vX,spm_vec(X{i}));
    end
else
    vX = [];
end
