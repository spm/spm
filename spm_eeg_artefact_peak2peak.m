function res = spm_eeg_artefact_peak2peak(S)
% Plugin for spm_eeg_artefact doing artefact detection based on peak-to-peak amplitude
% S            - input structure
% fields of S:
%    S.D       - M/EEG object
%    S.chanind - vector of indices of channels that this plugin will look at.
%                         
%    Additional parameters can be defined specific for each plugin.
%
% Output:
% res - 
%    If no input is provided the plugin returns a cfg branch for itself.
%
%    If input is provided the plugin returns a matrix of size D.nchannels x D.ntrials  
%    with zeros for clean channel/trials and ones for artefacts.
%__________________________________________________________________________
% Copyright (C) 2008-2017 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_artefact_peak2peak.m 7132 2017-07-10 16:22:58Z guillaume $


%-This part if for creating a config branch that plugs into spm_cfg_eeg_artefact
% Any parameters can be specified and they are then passed to the plugin
% when it's called.
%--------------------------------------------------------------------------
if nargin == 0
    threshold         = cfg_entry;
    threshold.tag     = 'threshold';
    threshold.name    = 'Threshold';
    threshold.strtype = 'r';
    threshold.num     = [1 1];
    threshold.help    = {'Threshold value to apply to all channels.'};

    peak2peak      = cfg_branch;
    peak2peak.tag  = 'peak2peak';
    peak2peak.name = 'Peak to peak amplitude';
    peak2peak.val  = {threshold};
    peak2peak.help = {''};
    
    res = peak2peak;
    
    return
end

SVNrev = '$Rev: 7132 $';

%-Startup
%--------------------------------------------------------------------------
spm('sFnBanner', mfilename, SVNrev);
spm('FigName','M/EEG peak to peak artefact detection');

if isequal(S.mode, 'mark')
    error('Only reject mode is supported by this plug-in.');
end

D = spm_eeg_load(S.D);

chanind = S.chanind;
threshold = S.threshold;
res = zeros(D.nchannels, D.ntrials);

%-Artefact detection
%--------------------------------------------------------------------------

spm_progress_bar('Init', D.ntrials, 'Trials checked');
if D.ntrials > 100, Ibar = floor(linspace(1, D.ntrials,100));
else Ibar = [1:D.ntrials]; end

for i = 1:D.ntrials
    res(chanind, i) = (squeeze(max(D(chanind, :, i), [], 2) - min(D(chanind, :, i), [], 2)))>threshold;
    if any(Ibar == i), spm_progress_bar('Set', i); end
end

spm_progress_bar('Clear');

spm('FigName','M/EEG peak to peak artefact detection: done');
