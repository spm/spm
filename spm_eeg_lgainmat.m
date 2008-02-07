function [L] = spm_eeg_lgainmat(D,Is)
% loads (memory maps) a gain matrix
% FORMAT [L] = spm_eeg_lgainmat(D,Is)
% D    - Data strcucture
% Is   - indices of vertices
%
% L    - Lead-field or gain matrix L(:,Is)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_eeg_lgainmat.m 1143 2008-02-07 19:33:33Z spm $

% get matrix
%--------------------------------------------------------------------------
fname = D.inv{D.val}.forward.gainmat;
try
    G = load(fname);
catch
    [p f] = fileparts(fname);
    try
        fname = fullfile(D.path,f);
        G     = load(fname);
    catch
        fname = fullfile(pwd,f);
        G     = load(fname);
    end
end
Gname = fieldnames(G);

% condition units and subsample
%--------------------------------------------------------------------------
L     = sparse(getfield(G, Gname{1}));
L     = spm_cond_units(L);
try
    L = L(:,Is);
end