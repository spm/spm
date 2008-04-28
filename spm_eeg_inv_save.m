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
% $Id: spm_eeg_inv_save.m 1490 2008-04-28 11:16:29Z vladimir $

% save D
%--------------------------------------------------------------------------

D.save;
