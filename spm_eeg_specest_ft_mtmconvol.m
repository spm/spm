function res = spm_eeg_specest_ft_mtmconvol(S, data, time)
% Plugin for spm_eeg_tf implementing spectral estimation using Fieldtrip's freqanalysis_mtmconvol
% FORMAT res = spm_eeg_specest_ft_mtmconvol(S, data, time)
%
% S                     - input structure
% fields of S:
%    S.taper       - taper to use ('hanning', 'rectwin', 'dpss', 'sine' or
%                    other possible inputs of 'window'
%    S.freqres     - frequency resolutions (plus-minus for each frequency, can
%                    be a vector with a value per frequency)
%    S.frequencies - vector of frequencies
%    S.timeres     - time resolution in ms (length of the sliding time-window)
%    S.timestep    - time step (in ms) to slide the time-window by.
%
% Output:
%  res -
%   If no input is provided the plugin returns a cfg branch for itself
%
%   If input is provided:
%      res.fourier - the complex output of wavelet transform (in the case
%                    of single taper)
%      res.pow     - power (in case of multiple tapers, phase is not computed)
%      res.time    - time axis
%      res.freq    - frequency axis
%______________________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak based on the code contributed by Krish Singh
% $Id: spm_eeg_specest_ft_mtmconvol.m 3876 2010-05-07 18:51:03Z vladimir $


%-This part if for creating a config branch that plugs into spm_cfg_eeg_tf
% Any parameters can be specified and they are then passed to the plugin
% when it's called.
%--------------------------------------------------------------------------
if nargin == 0
    timeres = cfg_entry;
    timeres.tag = 'timeres';
    timeres.name = 'Time resolution';
    timeres.strtype = 'r';
    timeres.num = [1 1];
    timeres.val = {400};
    timeres.help = {'Length of the sliding time window (in ms)'};
    
    timestep = cfg_entry;
    timestep.tag = 'timestep';
    timestep.name = 'Time step';
    timestep.strtype = 'r';
    timestep.num = [1 1];
    timestep.val = {50};
    timestep.help = {'Step to slide the time window by (in ms)'};
    
    freqres = cfg_entry;
    freqres.tag = 'freqres';
    freqres.name = 'Frequency resolution';
    freqres.strtype = 'r';
    freqres.num = [1 Inf];
    freqres.val = {2};
    freqres.help = {'Frequency resolution.',...
        'Note: 1 Hz resolution means plus-minus 1 Hz, i.e. 2 Hz badwidth',...
        'Either a single value or a vector of the same length as frequencies can be input'};
    
    taper = cfg_menu;
    taper.tag = 'taper';
    taper.name = 'Taper';
    taper.help = {'Save taper as well as power'};
    taper.labels = {'Hanning', 'Rectangular', 'DPSS', 'Sine'};
    taper.values = {'hanning', 'rectwin', 'dpss', 'sine'};
    taper.val = {'sine'};
    
    ft_mtmconvol = cfg_branch;
    ft_mtmconvol.tag = 'ft_mtmconvol';
    ft_mtmconvol.name = 'Fieldtrip multi-taper';
    ft_mtmconvol.val = {taper, timeres, timestep, freqres};
    
    res = ft_mtmconvol;
    
    return
elseif nargin < 3
    error('Three input arguments are required');
end

%-Defaults
%--------------------------------------------------------------------------
if ~isfield(S, 'taper')
    S.taper = 'dpss';
end

if ~isfield(S, 'timeres')
    S.timeres = 400;
end

timeres = 1e-3*S.timeres;

if ~isfield(S, 'timestep')
    S.timestep = 50;
end
timestep = 1e-3*S.timestep;

dt = time(end) - time(1);

if ~isfield(S, 'frequencies') || isempty(S.frequencies)
    S.frequencies = (1/dt):max(1/dt, floor(dt)/dt):48;
end

if ~isfield(S, 'freqres')
    S.freqres = max(1/timeres, floor(timeres)/timeres);
end

if length(S.freqres) == 1
    freqres = S.freqres*ones(1, length(S.frequencies));
elseif length(S.freqres) == length(S.frequencies)
    freqres = S.freqres;
else
    error('Frequency resolution should be either a scalar or a vector the same length as the number of frequencies.')
end

%-Data dimensions
%--------------------------------------------------------------------------
Nchannels = size(data, 1);
Nsamples = size(data, 2);
Nfrequencies = length(S.frequencies);

fsample = 1./diff(time(1:2));

%-Do the spectral analysis
%--------------------------------------------------------------------------
% This is a temporary solution to get the functionality that will soon
% be replaced by a low level specest function

ftraw = [];
ftraw.trial = {data};
ftraw.time = {time};
ftraw.fsample = fsample;

for i = 1:Nchannels
    ftraw.label{i, 1} = ['Ch' num2str(i)];
end

cfg = [];
cfg.output ='fourier';
cfg.taper = S.taper;
cfg.method          = 'mtmconvol';

% This sets the centers of frequency bins at the optimal locations based on
% the time window.
cfg.foi             = S.frequencies; % Frequency axis
numfoi              = length(cfg.foi);

% This means that the time resolution is the same for all frequencies
cfg.t_ftimwin       = zeros(1,numfoi);
cfg.t_ftimwin(:)    =  timeres; % Time resolution

cfg.tapsmofrq   = freqres;

cfg.toi=(time(1)+(timeres/2)):timestep:(time(end)-(timeres/2)-1/fsample); % Time axis

if ismember(cfg.taper, {'dpss', 'sine'}) && ~(all(cfg.tapsmofrq==cfg.tapsmofrq(1)) && all(cfg.t_ftimwin==cfg.t_ftimwin(1)))
    cfg.output ='pow';
end

freq = ft_freqanalysis(cfg, ftraw);

res = [];
res.freq = freq.freq;
res.time = freq.time;

if isfield(freq, 'fourierspctrm')
    if ndims(freq.fourierspctrm) == 4 && size(freq.fourierspctrm, 1)>1
        res.pow = spm_squeeze(mean(freq.fourierspctrm.*conj(freq.fourierspctrm), 1), 1);
    elseif ndims(freq.fourierspctrm) == 4
        res.fourier = spm_squeeze(freq.fourierspctrm, 1);
    else
        res.fourier = freq.fourierspctrm;
    end
else
    res.pow = freq.powspctrm;
end