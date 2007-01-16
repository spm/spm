function spm_eeg_inv_checkforward(D)
% checks forward model
% FORMAT spm_eeg_inv_checkforward(D)

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