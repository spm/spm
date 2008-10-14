function [L, D] = spm_eeg_lgainmat(D,Is)
% loads (memory maps) a gain matrix
% FORMAT [L] = spm_eeg_lgainmat(D,Is)
% D    - Data strcucture
% Is   - indices of vertices
%
% L    - Lead-field or gain matrix L(:,Is)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_eeg_lgainmat.m 2339 2008-10-14 18:39:21Z vladimir $

% get matrix
%--------------------------------------------------------------------------
try
    fname = D.inv{D.val}.forward.gainmat;
    try
        G = load(fname); % Absolute path
    catch
        try
            G = load(fullfile(D.path, fname)); % Relative path
        catch
            [p f x] = fileparts(fname);
            try
                fname = fullfile(D.path, [f x]); % Try in the same directory where D is
                G     = load(fname);
            catch
                try
                    fname = fullfile(pwd, [f x]); % Try in the current directory
                    G     = load(fname);
                end
            end
        end
    end
    modality = spm_eeg_modality_ui(D, 1);
    chanind = strmatch(modality, D.chantype, 'exact');
    chanind = setdiff(chanind, D.badchannels);
    Gname = fieldnames(G);
    L  = sparse(getfield(G, Gname{1}));
    if length(chanind)~=size(L, 1)
        error('Improper gain matrix');
    end
catch
    % That might be an overkill, but lets try
    D = spm_eeg_inv_forward(D, D.val);
    fname = D.inv{D.val}.forward.gainmat;
    G = load(fullfile(D.path, fname)); 
    Gname = fieldnames(G);
    L     = sparse(getfield(G, Gname{1}));
end

L     = spm_cond_units(L);

if nargin>1
    L = L(:,Is);
end