function res = spm_eeg_artefact_threshchan(S)
% Plugin for spm_eeg_artefact doing artefact detection by chanel thresholding.
% S                     - input structure
% fields of S:
%    S.D                - M/EEG object
%    S.chanind          - vector of indices of channels that this plugin will look at.
%                         
%    Additional parameters can be defined specific for each plugin
% Output:
%  res - 
%   If no input is provided the plugin returns a cfg branch for itself
%
%   If input is provided the plugin returns a matrix of size D.nchannels x D.ntrials  
%   with zeros for clean channel/trials and ones for artefacts.
%______________________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_artefact_threshchan.m 3262 2009-07-09 12:10:53Z vladimir $


%-This part if for creating a config branch that plugs into spm_cfg_eeg_artefact
% Any parameters can be specified and they are then passed to the plugin
% when it's called.
%--------------------------------------------------------------------------
if nargin == 0
    threshold = cfg_entry;
    threshold.tag = 'threshold';
    threshold.name = 'Threshold';
    threshold.strtype = 'r';
    threshold.num = [1 1];
    threshold.help = {'Threshold value to apply to all channels'};

    threshchan = cfg_branch;
    threshchan.tag = 'threshchan';
    threshchan.name = 'Threshold channels';
    threshchan.val = {threshold};
    
    res = threshchan;
    
    return
end

SVNrev = '$Rev: 3262 $';

%-Startup
%--------------------------------------------------------------------------
spm('sFnBanner', mfilename, SVNrev);
spm('FigName','M/EEG threshold channels');

D = spm_eeg_load(S.D);

chanind  = S.chanind;
threshold = S.threshold;
res = zeros(D.nchannels, D.ntrials);

%-Artefact detection
%--------------------------------------------------------------------------

spm_progress_bar('Init', D.ntrials, 'Trials checked');
if D.ntrials > 100, Ibar = floor(linspace(1, D.ntrials,100));
else Ibar = [1:D.ntrials]; end

for i = 1:D.ntrials
    res(chanind, i) = squeeze(max(abs(D(chanind, :, i)), [], 2))>threshold;
    if ismember(i, Ibar), spm_progress_bar('Set', i); end
end

spm_progress_bar('Clear');

spm('FigName', 'M/EEG threshold channels: done');