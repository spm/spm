function Dtf = spm_eeg_ft_multitaper_powermap(S)
% Function for computing spectra of stationary data using multitaper
%
% FORMAT D = spm_eeg_ft_multitaper_powermap(S)
%
% S           - input structure (optional)
% (optional) fields of S:
%   S.D       - filename, or M/EEG object
%   S.freqwin  - frequency window (Hz)
%   S.freqres  - frequency resolution (Hz)
%   S.log      - compute log of power (1 - yes, 0- no)
%
% Output
%   Dtf - dataset with the spectra where the frequency dimension has been
%         converted to time dimension as a hack to allow export to images
% Disclaimer: this code is provided as an example and is not guaranteed to work
% with data on which it was not tested. If it does not work for you, feel
% free to improve it and contribute your improvements to the MEEGtools toolbox
% in SPM (http://www.fil.ion.ucl.ac.uk/spm)
%
% _______________________________________________________________________
% Copyright (C) 2009 Institute of Neurology, UCL

% Vladimir Litvak
% $Id: spm_eeg_ft_multitaper_powermap.m 5986 2014-05-15 09:36:55Z vladimir $
 
%%
SVNrev = '$Rev: 5986 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','Fieldtrip multitaper power map'); 

%%
%-Test for the presence of required Matlab toolbox
%--------------------------------------------------------------------------
if ~license('test','signal_toolbox')
    error('This function requires the Signal Processing Toolbox.');
end

%% ============ Load SPM EEG file and verify consistency
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

%-Configure the spectral analysis
%--------------------------------------------------------------------------
if ~isfield(S, 'timewin')
    S.timewin = spm_input('Time window (ms)', '+1', 'r', num2str(1000*[D.time(1) D.time(end)]), 2);
end

if ~isfield(S, 'freqwin')
    S.freqwin = spm_input('Frequency window (Hz)', '+1', 'r', '0 90', 2);
end

if ~isfield(S, 'freqres')
    S.freqres = spm_input('Frequency resolution (Hz)', '+1', 'r', '1', 1);
end

if ~isfield(S, 'log')
    S.log = spm_input('Compute log(power)', '+1', 'yes|no', [1 0], 1);
end

data = D.ftraw;
%%

%-Run the Fieldtrip code
%--------------------------------------------------------------------------
cfg = [];
cfg.latency = 1e-3*S.timewin;
cfg.keeptrials = 'yes';

data = ft_timelockanalysis(cfg, data);


cfg = [];
cfg.output ='pow';
cfg.keeptrials = 'yes';
cfg.taper = 'dpss';
cfg.channel = D.chanlabels(D.indchantype('MEEG', 'GOOD'));
cfg.method          = 'mtmfft';
cfg.foilim          = S.freqwin; 
cfg.tapsmofrq       = S.freqres;

freq = ft_freqanalysis(cfg, data);


%-Create a dummy Fieldtrip structure pretending frequency is time
%--------------------------------------------------------------------------

dummy = [];
dummy.dimord = 'rpt_chan_time';
dummy.time   = 1e-3*freq.freq;

if length(dummy.time)>1
    dummy.fsample = 1/(dummy.time(2)-dummy.time(1));
else
    dummy.fsample = 1;
end

dummy.label = freq.label;
dummy.trial  = freq.powspctrm;

megind = spm_match_str(freq.label, D.chanlabels(D.indchantype('MEG', 'GOOD')));
if ~isempty(megind)
    dummy.trial(:, megind, :) = 1e30*dummy.trial(:, megind, :);
end

if S.log
    dummy.trial = log(dummy.trial);
end

%-Save the result in SPM dataset
%--------------------------------------------------------------------------
Dtf  = spm_eeg_ft2spm(dummy, fullfile(D.path, ['TF' D.fname])); 

%-Copy some additional information from the original file
%--------------------------------------------------------------------------
Dtf  = conditions (Dtf, ':', D.conditions);

[sel1, sel2] = spm_match_str(Dtf.chanlabels, D.chanlabels);

Dtf = chantype(Dtf, sel1, chantype(D, sel2));
Dtf = badchannels(Dtf, sel1, badchannels(D, sel2));
Dtf = coor2D(Dtf, ':', coor2D(D, ':'));

Dtf = badtrials(Dtf, badtrials(D), 1);
Dtf = history(Dtf, history(D));

%-Update history
%--------------------------------------------------------------------------
Dtf = history(Dtf, mfilename, S);

%-Save
%--------------------------------------------------------------------------
save(Dtf);



