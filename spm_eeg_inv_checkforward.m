function spm_eeg_inv_checkforward(D)
% checks forward model
% FORMAT spm_eeg_inv_checkforward(D)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jeremie Mattout
% $Id: spm_eeg_inv_checkforward.m 1437 2008-04-17 10:34:39Z christophe $

% SPM data structure
%==========================================================================
try
    forward = D.inv{D.ival}.forward;
    disp(forward)
catch
    warndlg('please create forward model')
    return
end

% show electrodes
%--------------------------------------------------------------------------
if strcmp(D.inv{D.ival}.method,'ECD')
    
    % SPM graphics figure
    %----------------------------------------------------------------------
    Fgraph = spm_figure('GetWin','Graphics');
    spm_figure('Clear',Fgraph)

    spm_eeg_inv_displScEl(forward.head(end),forward.electrodes);
end
