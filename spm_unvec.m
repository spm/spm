function [X] = spm_vec(vX,X)
% unvectorises a vectorised array 
% FORMAT [X] = spm_vec(vX,X);
% X  - numeric, cell or stucture array
% vX - spm_vec(X)
%
% i.e. X = spm_unvec(spm_vec(X),X)
%
% see spm_vec
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id: spm_vec.m 184 2005-05-31 13:23:32Z karl $

% fill in structure arrays
%--------------------------------------------------------------------------
if isstruct(X)
    f = fieldnames(X);
    for i = 1:length(f)
        c          = {X.(f{i})};
        n          = 1:length(spm_vec(c));
        c          = spm_unvec(vX(n),c);
        [X.(f{i})] = deal(c{:});
        vX(n)      = [];
    end
    return
end

% fill in cells arrays
%--------------------------------------------------------------------------
if iscell(X)
    for i = 1:length(X(:))
        n     = 1:length(spm_vec(X{i}));
        X{i}  = spm_unvec(vX(n),X{i});
        vX(n) = [];
    end
    return
end

% reshape numerical arrays
%--------------------------------------------------------------------------
if isnumeric(X)
    X(:)  = vX;
else
    X     = [];
end
