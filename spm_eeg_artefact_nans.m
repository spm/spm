function res = spm_eeg_artefact_nans(S)
% Plugin for spm_eeg_artefact doing NaN detection.
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
% Copyright (C) 2011 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_artefact_nans.m 4350 2011-06-13 16:31:02Z vladimir $


%-This part if for creating a config branch that plugs into spm_cfg_eeg_artefact
% Any parameters can be specified and they are then passed to the plugin
% when it's called.
%--------------------------------------------------------------------------
if nargin == 0    
    nans = cfg_branch;
    nans.tag = 'nans';
    nans.name = 'Detect NaNs';
    nans.val = {};
    
    res = nans;
    
    return
end

SVNrev = '$Rev: 4350 $';

%-Startup
%--------------------------------------------------------------------------
spm('sFnBanner', mfilename, SVNrev);
spm('FigName','M/EEG NaN detection');

D = spm_eeg_load(S.D);

chanind  =  S.chanind;
res = zeros(D.nchannels, D.ntrials);

%-Artefact detection
%--------------------------------------------------------------------------

spm_progress_bar('Init', D.ntrials, 'Trials checked');
if D.ntrials > 100, Ibar = floor(linspace(1, D.ntrials,100));
else Ibar = [1:D.ntrials]; end

for i = 1:D.ntrials
    for j = 1:length(chanind)       
        if any(isnan(squeeze(D(chanind(j), :, i))))
            res(chanind(j), i) = 1;
        end
    end
    if ismember(i, Ibar), spm_progress_bar('Set', i); end
end

spm_progress_bar('Clear');

spm('FigName','M/EEG NaN detection: done');