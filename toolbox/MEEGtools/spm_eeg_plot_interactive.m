% Demo script for interactive plotting in Fieldtrip
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_plot_interactive.m 1673 2008-05-16 15:32:22Z vladimir $

D = spm_eeg_load;

data = D.ftraw;

if D.ntrials > 1
    clb = D.conditions;
    ind = spm_input('Select trial',1, 'm', sprintf('%s|', clb{:}),1:D.ntrials);
    
    data.trial = data.trial(ind);
    data.time =  data.time(ind);
end

data = ft_checkdata(data, 'datatype', 'timelock');

cfg = [];
cfg.interactive = 'yes';

switch spm_eeg_modality_ui(D)
    case 'EEG'
        chanind = strmatch('EEG', D.chantype);
        cfg.elec = D.sensors('EEG');
        data.elec = cfg.elec;
    case 'MEG'
        chanind = strmatch('MEG', D.chantype);
        cfg.grad = D.sensors('MEG');
        data.grad = cfg.grad;
end

cfg.channel = data.label(chanind);
%%
figure;
ft_multiplotER(cfg, data);



