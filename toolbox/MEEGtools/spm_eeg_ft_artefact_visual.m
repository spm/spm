function D = spm_eeg_ft_artefact_visual(S)
% Function for interactive artefact rejection using Fieldtrip
%
% Disclaimer: this code is provided as an example and is not guaranteed to work
% with data on which it was not tested. If it does not work for you, feel
% free to improve it and contribute your improvements to the MEEGtools toolbox
% in SPM (http://www.fil.ion.ucl.ac.uk/spm)
%
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_ft_artefact_visual.m 5674 2013-10-09 10:00:26Z vladimir $

[Finter,Fgraph,CmdLine] = spm('FnUIsetup', 'Fieldtrip visual artefact rejection',0);

if nargin == 0
    S = [];
end

try
    D = S.D;
catch
    D = spm_select(1, '\.mat$', 'Select EEG mat file');
    S.D = D;
end

if ischar(D)
    try
        D = spm_eeg_load(D);
    catch
        error(sprintf('Trouble reading file %s', D));
    end
end

if isequal(D.type, 'continuous')
    error('Artefact detection can only be applied to epoched data');
end

trlind = find(~badtrials(D, ':'));

data = ftraw(D, ':', ':', trlind);

data.trialinfo = [1:length(trlind)]';

cfg=[];

if isfield(S, 'method')
    cfg.method = S.method;
else
    cfg.method =  spm_input('What method?','+1', 'm', 'summary|channel|trial', strvcat('summary', 'channel', 'trial'));
end

if isfield(S, 'latency')
    cfg.latency = S.latency;
else
    cfg.latency = 1e-3*spm_input('PST ([start end] in ms):', '+1', 'r', num2str(1e3*[data.time{1}(1) data.time{1}(end)]), 2);
end

cfg.keepchannel = 'no';

if isfield(S, 'channel')
    cfg.channel = S.channel;
else
    [junk, chanind] = spm_eeg_modality_ui(D, 0, 1);
    cfg.channel = data.label(chanind);
end

data = ft_rejectvisual(cfg, data);

% Figure out based on the output of FT function what trials and channels to
% reject
trlsel = ones(1, length(trlind));
trlsel(data.trialinfo (:, 1)) = 0;

D = badtrials(D, trlind, trlsel);

badchan = setdiff(cfg.channel, data.label);
if ~isempty(badchan)
    badchanind = spm_match_str(D.chanlabels, badchan);
    D = badchannels(D, badchanind, 1);
end

if ~isfield(S, 'save') || S.save
    save(D);
end