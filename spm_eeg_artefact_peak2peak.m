function res = spm_eeg_artefact_peak2peak(S)
% Plugin for spm_eeg_artefact doing artefact detection based on peak-to-peak amplitude.
%______________________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_artefact_peak2peak.m 3258 2009-07-08 17:46:54Z vladimir $


if nargin == 0
    threshold = cfg_entry;
    threshold.tag = 'thresholdval';
    threshold.name = 'Threshold';
    threshold.strtype = 'r';
    threshold.num = [1 1];
    threshold.help = {'Threshold value to apply to all channels'};

    peak2peak = cfg_branch;
    peak2peak.tag = 'peak2peak';
    peak2peak.name = 'Peak to peak amplitude';
    peak2peak.val = {threshold};
    
    res = peak2peak;
    
    return
end

SVNrev = '$Rev: 3258 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG peak to peak artefact detection');

D = spm_eeg_load(S.D);

channels  = S.channels;
threshold = S.threshold;
res = zeros(D.nchannels, D.ntrials);

%-Artefact detection
%--------------------------------------------------------------------------

spm_progress_bar('Init', D.ntrials, 'Trials thresholded');
if D.ntrials > 100, Ibar = floor(linspace(1, D.ntrials,100));
else Ibar = [1:D.ntrials]; end

for i = 1:D.ntrials
    res(channels, i) = (squeeze(max(D(channels, :, i), [], 2) - min(D(channels, :, i), [], 2)))>threshold;
    if ismember(i, Ibar), spm_progress_bar('Set', i); end
end

spm_progress_bar('Clear');

spm('FigName','M/EEG peak to peak artefact detection: done');