function spm_eeg_inv_checkforward(D)
% checks forward model
% FORMAT spm_eeg_inv_checkforward(D)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jeremie Mattout
% $Id: spm_eeg_inv_checkforward.m 1488 2008-04-27 14:11:48Z vladimir $

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
    
    % SPM graphics figure
    %----------------------------------------------------------------------
    Fgraph = spm_figure('GetWin','Graphics');
    spm_figure('Clear',Fgraph)

    spm_eeg_inv_displScEl(forward.head(end),forward.electrodes);
end
