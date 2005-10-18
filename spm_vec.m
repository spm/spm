function [vX] = spm_vec(X)
% vectorises a numeric, cell or structure array
% FORMAT [vX] = spm_vec(X);
% X  - numeric, cell or stucture array
% vX - vec(X)
%__________________________________________________________________________
%
% e.g.:
% spm_vec({eye(2) 3}) = [1 0 0 1 3]'
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience
 
% Karl Friston
% $Id: spm_vec.m 258 2005-10-18 18:21:07Z karl $


% initialise vX
%--------------------------------------------------------------------------
vX = [];

% vectorise structure into cell arrays
%--------------------------------------------------------------------------
if isstruct(X)
    f = fieldnames(X);
    X = X(:);
    for i = 1:length(f)
            vX = [vX; spm_vec({X.(f{i})})];
    end
    return
end
 
% vectorise cells into numerical arrays
%--------------------------------------------------------------------------
if iscell(X)
    X     = X(:);
    for i = 1:length(X)
         vX = [vX; spm_vec(X{i})];
    end
    return
end
 
% vectorise numerical arrays
%--------------------------------------------------------------------------
if isnumeric(X)
    vX = X(:);
end

