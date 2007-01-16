function [D] = spm_eeg_inv_results_ui(varargin)
% GUI for contrast of evoked responses and power for an MEG-EEG model
% FORMAT [D] = spm_eeg_inv_results_ui(D,val)
% Sets:
%
%     D.contrast.woi   - time (ms) window of interest
%     D.contrast.fboi  - freq (Hz) window of interest
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Jeremie Mattout
% $Id: spm_eeg_inv_results_ui $

% initialise
%--------------------------------------------------------------------------
[D,val] = spm_eeg_inv_check(varargin{:});
try
    D.inv{val}.inverse.J;
catch
    warndlg('Please invert this model first');
    return
end

% get time window
%--------------------------------------------------------------------------
woi    = spm_input('Time window (ms)','+1','r',[100 200]);
woi    = sort(woi);
D.inv{val}.contrast.woi = round([woi(1) woi(end)]);

% get frequency window
%--------------------------------------------------------------------------
fboi  = spm_input('Frequency [band] of interest (Hz)','+1','r',0);
fboi  = sort(fboi);
D.inv{val}.contrast.fboi = round([fboi(1) fboi(end)]);

% evaluate contrast
%--------------------------------------------------------------------------
D     = spm_eeg_inv_results(D);
