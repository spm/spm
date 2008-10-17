function [L,D] = spm_eeg_lgainmat(D,Is)
% loads (memory maps) a gain matrix
% FORMAT [L,D] = spm_eeg_lgainmat(D,Is)
% D    - Data structure
% Is   - indices of vertices
%
% L    - Lead-field or gain matrix L(:,Is)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_eeg_lgainmat.m 2352 2008-10-17 11:53:53Z karl $
 
% get gain or lead-feild matrix
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
                fname = fullfile(D.path, [f x]); % Try D's directory
                G     = load(fname);
            catch
                try
                    fname = fullfile(pwd, [f x]); % Try current directory
                    G     = load(fname);
                end
            end
        end
    end
    modality = spm_eeg_modality_ui(D,1);
    chanind  = strmatch(modality, D.chantype, 'exact');
    chanind  = setdiff(chanind, D.badchannels);
    Gname    = fieldnames(G);
    L        = sparse(getfield(G,Gname{1}));
    if length(chanind) ~= size(L,1)
        error('Gain matrix has an incorrect number of channels');
    end
catch
 
    % create a new lead-field matrix
    %----------------------------------------------------------------------
    D     = spm_eeg_inv_forward(D, D.val);
    fname = D.inv{D.val}.forward.gainmat;
    G     = load(fullfile(D.path, fname));
    Gname = fieldnames(G);
    L     = sparse(getfield(G, Gname{1}));
end
 
 
% condition the scaling of the lead-field
%--------------------------------------------------------------------------
L     = spm_cond_units(L);
 
% retain selected sources if necessary
%--------------------------------------------------------------------------
if nargin > 1
    L = L(:,Is);
end
