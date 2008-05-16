% Demo script for interactive artefact rejection using Fieldtrip
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_ft_artefact_visual.m 1673 2008-05-16 15:32:22Z vladimir $

D = spm_eeg_load;

data = D.ftraw(0);

cfg=[];
cfg.method =  spm_input('What method?','+1', 'm', 'summary|channel|trial', strvcat('summary', 'channel', 'trial'));
cfg.keepchannel = 'yes';


switch spm_eeg_modality_ui(D)
    case 'EEG'
        chanind = strmatch('EEG', D.chantype);
    case 'MEG'
        chanind = strmatch('MEG', D.chantype);
end

cfg.channel = data.label(chanind);

[data, chansel, trlsel] = ft_rejectvisual(cfg, data);

D = reject(D, 1:D.ntrials, ~trlsel);
D = badchannels(D, 1:D.nchannels, ~chansel);

save(D);