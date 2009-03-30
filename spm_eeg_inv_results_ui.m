function [D] = spm_eeg_inv_results_ui(varargin)
% GUI for contrast of evoked responses and power for an MEG-EEG model
% FORMAT [D] = spm_eeg_inv_results_ui(D,val)
% Sets:
%
%     D.contrast.woi   - time (ms) window of interest
%     D.contrast.fboi  - freq (Hz) window of interest
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_eeg_inv_results_ui.m 3000 2009-03-30 16:15:00Z will $

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
D.inv{val}.contrast.woi = round([min(woi) max(woi)]);

% get frequency window
%--------------------------------------------------------------------------
fboi  = spm_input('Frequency [band] of interest (Hz)','+1','r',0);
D.inv{val}.contrast.fboi = round([min(fboi) max(fboi)]);

% induced or evoked? (not for fusion)
%--------------------------------------------------------------------------
if strcmp(D.type, 'single')
    str  = {'evoked','induced'};
    type = spm_input('Power','+1','b',str,[],1);
else
    type = 'evoked';
end
D.inv{val}.contrast.type = type;

% evaluate contrast
%--------------------------------------------------------------------------
D     = spm_eeg_inv_results(D);
