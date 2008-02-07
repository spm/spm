function spm_eeg_inv_save(D)
% saves EEG/MEG D.mat data structure in fullfile(D.path,D.fname)
% FORMAT spm_eeg_inv_save(S)
%
% D  - D.mat data structure
%__________________________________________________________________________
%
% Ensures matlab version backwards compatibility
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_eeg_inv_save.m 1143 2008-02-07 19:33:33Z spm $

% save D
%--------------------------------------------------------------------------
if spm_matlab_version_chk('7.1') >= 0
    save(fullfile(D.path, D.fname), '-V6', 'D');
else
    save(fullfile(D.path, D.fname), 'D');
end

