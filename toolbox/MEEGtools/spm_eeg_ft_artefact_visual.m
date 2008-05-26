% Demo script for interactive artefact rejection using Fieldtrip
%
% Disclaimer: this code is provided as an example and is not guaranteed to work
% with data on which it was not tested. If it does not work for you, feel
% free to improve it and contribute your improvements to the MEEGtools toolbox
% in SPM (http://www.fil.ion.ucl.ac.uk/spm)
%
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_ft_artefact_visual.m 1726 2008-05-26 16:45:55Z vladimir $

D = spm_eeg_load;

data = D.ftraw(0);

trlind = 1:length(data.trial);
data.cfg.trl(:, 4) = trlind;

cfg=[];
cfg.method =  spm_input('What method?','+1', 'm', 'summary|channel|trial', strvcat('summary', 'channel', 'trial'));
cfg.keepchannel = 'no';


switch spm_eeg_modality_ui(D)
    case 'EEG'
        chanind = strmatch('EEG', D.chantype);
    case 'MEG'
        chanind = strmatch('MEG', D.chantype);
end

cfg.channel = data.label(chanind);

data = ft_rejectvisual(cfg, data);

% Figure out based on the output of FT function what trials and channels to
% reject
trlsel(trlind) = 0;
trlsel(data.cfg.trl(:, 4)) = 1;
D = reject(D, 1:D.ntrials, ~trlsel);

badchan = setdiff(cfg.channel, data.label);
if ~isempty(badchan)
    badchanind = spm_match_str(D.chanlabels, badchan);
    D = badchannels(D, badchanind, 1);
end

save(D);