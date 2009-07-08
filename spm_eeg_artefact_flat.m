function res = spm_eeg_artefact_flat(S)
% Plugin for spm_eeg_artefact doing flat channel detection.
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
% $Id: spm_eeg_artefact_flat.m 3258 2009-07-08 17:46:54Z vladimir $


%-This part if for creating a config branch that plugs into spm_cfg_eeg_artefact
% Any parameters can be specified and they are then passed to the plugin
% when it's called.
%--------------------------------------------------------------------------
if nargin == 0
    threshold = cfg_entry;
    threshold.tag = 'threshold';
    threshold.name = 'Threshold';
    threshold.strtype = 'r';
    threshold.val = {0};
    threshold.num = [1 1];
    threshold.help = {'Threshold for difference between adjacent samples'};

    
    seqlength = cfg_entry;
    seqlength.tag = 'seqlength';
    seqlength.name = 'Flat segment length';
    seqlength.strtype = 'r';
    seqlength.num = [1 1];
    seqlength.val = {4};
    seqlength.help = {'Minimal number of adjacent samples with the same value to reject.'};

    flat = cfg_branch;
    flat.tag = 'flat';
    flat.name = 'Flat segments';
    flat.val = {threshold, seqlength};
    
    res = flat;
    
    return
end

SVNrev = '$Rev: 3258 $';

%-Startup
%--------------------------------------------------------------------------
spm('sFnBanner', mfilename, SVNrev);
spm('FigName','M/EEG flat data detection');

D = spm_eeg_load(S.D);

chanind  =  S.chanind;
threshold = S.threshold;
seqlength = S.seqlength;
res = zeros(D.nchannels, D.ntrials);

%-Artefact detection
%--------------------------------------------------------------------------

spm_progress_bar('Init', D.ntrials, 'Trials checked');
if D.ntrials > 100, Ibar = floor(linspace(1, D.ntrials,100));
else Ibar = [1:D.ntrials]; end

for i = 1:D.ntrials
    for j = 1:length(chanind)
        tmp  = find(abs(diff(squeeze(D(chanind(j), :, i)), [], 2))>threshold);
        if max(diff([0 tmp D.nsamples]))>=seqlength
            res(chanind(j), i) = 1;
        end
    end
    if ismember(i, Ibar), spm_progress_bar('Set', i); end
end

spm_progress_bar('Clear');

spm('FigName','M/EEG flat data detection: done');