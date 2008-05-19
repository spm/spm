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
% $Id: spm_eeg_ft_artefact_visual.m 1680 2008-05-19 11:18:19Z vladimir $

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