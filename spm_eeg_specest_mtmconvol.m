function res = spm_eeg_specest_mtmconvol(S, data, time)
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
% Copyright (C) 2011-2013 Wellcome Trust Centre for Neuroimaging

% $Id: spm_eeg_specest_mtmconvol.m 6612 2015-11-27 18:43:15Z vladimir $


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
    
    mtmconvol = cfg_branch;
    mtmconvol.tag = 'mtmconvol';
    mtmconvol.name = 'Multi-taper';
    mtmconvol.val = {taper, timeres, timestep, freqres};
    
    res = mtmconvol;
    
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

dt = time(end) - time(1) + diff(time(1:2));


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
fsample = 1./diff(time(1:2));

df = unique(diff(S.frequencies));
if length(df) == 1  
    pad = ceil(dt*df)/df;
else
    pad = [];
end

% Correct the time step to the closest multiple of the sampling interval to
% keep the time axis uniform
timestep = round(fsample*timestep)/fsample;

timeoi=(time(1)+(timeres/2)):timestep:(time(end)-(timeres/2)-1/fsample); % Time axis

[spectrum,ntaper,freqoi,timeoi] = ft_specest_mtmconvol(data, time, 'taper', S.taper, 'timeoi', timeoi, 'freqoi', S.frequencies,...
    'timwin', repmat(timeres, 1, length(S.frequencies)), 'tapsmofrq', freqres, 'pad', pad, 'verbose', 0);

% To prevent differences between files due to numerical imprecision
if length(freqoi)==length(S.frequencies) && max(abs(freqoi-S.frequencies))<0.1
    freqoi = S.frequencies;
end

% This in principle should not happen
if length(unique(diff(timeoi)))>1 && max(abs(diff(unique(diff(timeoi)))))>1e-6
    error('Non-uniform output time axis')
end

if length(unique(diff(freqoi)))>1 && max(abs(diff(unique(diff(freqoi)))))>1e-2
    error('Non-uniform output frequency axis')
end

res = [];
res.freq = freqoi;
res.time = timeoi;

if ndims(spectrum) == 4 && size(spectrum, 1)>1
    res.pow = spm_squeeze(nanmean(spectrum.*conj(spectrum), 1), 1);
elseif ndims(spectrum) == 4
    res.fourier = spm_squeeze(spectrum, 1);
else
    res.fourier = spectrum;
end
end

function y = nansum(x, dim)

x(isnan(x)) = 0;
if nargin==1
  y = sum(x);
else
  y = sum(x,dim);
end
end

function y = nanmean(x, dim)
N = sum(~isnan(x), dim);
y = nansum(x, dim) ./ N;
end