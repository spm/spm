function Dtf = spm_eeg_ft_multitaper_tf(S)
% Function for computing time-frequency decomposition using multitaper
%
% WARNING: This function uses some quite specific settings and is not generic. It is
% just an example of how Fieldtrip spectral analysis can be combined with
% SPM
%
% FORMAT D = spm_eeg_ft_multitaper_tf(S)
%
% S           - input structure (optional)
% (optional) fields of S:
%   S.D       - filename, or M/EEG object
%   S.channels - list of channel labels
%   S.pretrig  - time to start TF analysis in PST (ms)
%   S.posttrig - time to end TF analysis in PST (ms)
%   S.timewin  - time window (resolution) in ms
%   S.timestep - time step in ms
%   S.freqwin  - frequency window (Hz)
%
% Disclaimer: this code is provided as an example and is not guaranteed to work
% with data on which it was not tested. If it does not work for you, feel
% free to improve it and contribute your improvements to the MEEGtools toolbox
% in SPM (http://www.fil.ion.ucl.ac.uk/spm)
%
% _______________________________________________________________________
% Copyright (C) 2008 Institute of Neurology, UCL

% Vladimir Litvak
% $Id: spm_eeg_ft_multitaper_tf.m 3623 2009-12-09 09:40:36Z vladimir $
 
%%
SVNrev = '$Rev: 3623 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','Fieldtrip multitaper TF'); 

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
if ~isfield(S, 'channels')
    switch spm_input('What channels','+1','all|meeg|select', strvcat('all', 'meeg', 'select'), 'all')
        case 'all'
            S.channels = D.chanlabels;
        case 'meeg'
            S.channels = D.chanlabels(D.meegchannels);
        case 'select'
            [selection, ok]= listdlg('ListString', D.chanlabels, 'SelectionMode', 'multiple' ,'Name', 'Select channels' , 'ListSize', [400 300]);
            if ~ok
                return;
            end
            S.channels = D.chanlabels(selection);
    end
end

if ~isfield(S, 'pretrig')
    S.pretrig = spm_input('Start in PST [ms]', '+1', 'r', '', 1);
end

if ~isfield(S, 'posttrig')
    S.posttrig = spm_input('End in PST [ms]', '+1', 'r', '', 1);
end

if ~isfield(S, 'timewin')
    S.timewin = spm_input('Time window (ms)', '+1', 'r', '400', 1);
end

if ~isfield(S, 'timestep')
    S.timestep = spm_input('Time step (ms)', '+1', 'r', '50', 1);
end

if ~isfield(S, 'freqwin')
    S.freqwin = spm_input('Frequency window (Hz)', '+1', 'r', '0 90', 2);
end

data = D.ftraw;

prestim = 1e-3*S.pretrig;
poststim = 1e-3*S.posttrig;
%%
timewin = 1e-3*S.timewin;
step = 1e-3*S.timestep;

%-Run the Fieldtrip code
%--------------------------------------------------------------------------

cfg = [];
cfg.output ='pow';
cfg.keeptrials = 'yes';
cfg.taper = 'dpss';
cfg.channel = S.channels;
cfg.method          = 'mtmconvol';

% This sets the centers of frequency bins at the optimal locations based on
% the time window.
cfg.foi             = (1/timewin):(1/timewin):S.freqwin(2); % Frequency axis
cfg.foi             = cfg.foi(cfg.foi>=S.freqwin(1));
numfoi              = length(cfg.foi);

% This means that the time resolution is the same for all frequencies
cfg.t_ftimwin       = zeros(1,numfoi);
cfg.t_ftimwin(:)    = timewin; % Time resolution

% This part is about frequency resolution
cfg.tapsmofrq       = zeros(1,numfoi); % Initialize to zero
cfg.tapsmofrq(:)    = 1/timewin; % Set initial resolution to 1/timewin (i.e. 2.5 Hz) for all frequencis
% Here it sets the resolution for frequencies above 10*(1/timewin) (25 Hz)
% to 0.1 times the frequency. This means that at up to 25 Hz the resolution
% is fixed and then it starts slowly increasing up to 10Hz in each
% direction for 100 Hz. If you comment out this line, you'll have fixed
% time resolution for all frequencies.
cfg.tapsmofrq(cfg.foi>10*(1/timewin))    = 0.1*cfg.foi(cfg.foi>10*(1/timewin));
cfg.tapsmofrq(cfg.foi>50)                = 5;
% figure; plot(cfg.foi, cfg.tapsmofrq);xlabel('frequency (Hz)');ylabel('frequency resolution (Hz)')
% This is the time axis. The window with the width defined above (400 msec)
% is moved every time by 'step' (100 ms). The earliest you can start is
% half the time window from the start of the data. Otherwise your time
% window will overlap a segment with no data and you will get NaNs in the
% output. The same idea at the end.
cfg.toi=(prestim+(timewin/2)):step:(poststim-(timewin/2)-1/data.fsample); % Time axis

cfg.toi = cfg.toi((cfg.toi-(timewin/2)>=D.time(1)) & (cfg.toi+timewin/2+1/data.fsample)<=D.time(end));

freq = ft_freqanalysis(cfg, data);

%-Save the result in SPM dataset
%--------------------------------------------------------------------------
Dtf  = spm_eeg_ft2spm(freq, fullfile(D.path, ['TF' D.fname])); 

%-Copy some additional information from the original file
%--------------------------------------------------------------------------
Dtf  = events(Dtf, [], events(D, []));
Dtf  = conditions (Dtf, [], D.conditions);

[sel1, sel2] = spm_match_str(Dtf.chanlabels, D.chanlabels);

Dtf = chantype(Dtf, sel1, chantype(D, sel2));
Dtf = badchannels(Dtf, sel1, badchannels(D, sel2));
Dtf = coor2D(Dtf, [], coor2D(D, []));

Dtf = reject(Dtf, [], reject(D));
Dtf = history(Dtf, history(D));

%-Update history
%--------------------------------------------------------------------------
Dtf = history(Dtf, mfilename, S);

%-Save
%--------------------------------------------------------------------------
save(Dtf);



