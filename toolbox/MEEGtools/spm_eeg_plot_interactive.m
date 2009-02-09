% Demo script for interactive plotting in Fieldtrip
%
% Disclaimer: this code is provided as an example and is not guaranteed to work
% with data on which it was not tested. If it does not work for you, feel
% free to improve it and contribute your improvements to the MEEGtools toolbox
% in SPM (http://www.fil.ion.ucl.ac.uk/spm)
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_plot_interactive.m 2720 2009-02-09 19:50:46Z vladimir $

D = spm_eeg_load;

data = D.ftraw;

if D.ntrials > 1
    clb = D.conditions;
    ind = spm_input('Select trial',1, 'm', sprintf('%s|', clb{:}),1:D.ntrials);
    
    data.trial = data.trial(ind);
    data.time =  data.time(ind);
end

data = ft_timelockanalysis([], data);

cfg = [];
cfg.interactive = 'yes';

switch D.modality
    case 'EEG'
        chanind = strmatch('EEG', D.chantype);
        cfg.elec = D.sensors('EEG');
        cfg.rotate = 0;
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



