function spm_eeg_inv_checkforward(D)
% checks forward model
% FORMAT spm_eeg_inv_checkforward(D)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_eeg_inv_checkfoward.m 1131 2008-02-06 11:17:09Z spm $

% SPM data structure
%==========================================================================
try
    forward = D.inv{D.val}.forward;
    disp(forward)
catch
    warndlg('please create forward model')
    return
end

% show electrodes
%--------------------------------------------------------------------------
if strcmp(D.inv{D.val}.method,'ECD') 
    spm_eeg_inv_displScEl(forward.head(end),forward.electrodes);
end
