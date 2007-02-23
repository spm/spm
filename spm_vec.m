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
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience
 
% Karl Friston
% $Id: spm_vec.m 744 2007-02-23 13:51:48Z karl $

% initialise X and vX
%--------------------------------------------------------------------------
X     = varargin;
if length(X) == 1
    X = X{1};
end
vX    = [];

% vectorise structure into cell arrays
%--------------------------------------------------------------------------
if isstruct(X)
    f = fieldnames(X);
    X = X(:);
    for i = 1:length(f)
            vX = cat(1,vX,spm_vec({X.(f{i})}));
    end
    return
end
 
% vectorise cells into numerical arrays
%--------------------------------------------------------------------------
if iscell(X)
    X     = X(:);
    for i = 1:length(X)
         vX = cat(1,vX,spm_vec(X{i}));
    end
    return
end
 
% vectorise numerical arrays
%--------------------------------------------------------------------------
if isnumeric(X) || islogical(X)
    vX = X(:);
end

