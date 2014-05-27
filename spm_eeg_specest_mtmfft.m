function res = spm_eeg_specest_mtmfft(S, data, time)
% Plugin for spm_eeg_tf implementing spectral estimation using Fieldtrip's freqanalysis_mtmconvol
% FORMAT res = spm_eeg_specest_mtmfft(S, data, time)
%
% S                     - input structure
% fields of S:
%    S.taper       - taper to use ('hanning', 'rectwin', 'dpss', 'sine' or
%                    other possible inputs of 'window'
%    S.freqres     - frequency resolutions (plus-minus for each frequency, can
%                    be a vector with a value per frequency)
%    S.frequencies - vector of frequencies%
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
% Copyright (C) 2011-2014 Wellcome Trust Centre for Neuroimaging

% $Id: spm_eeg_specest_mtmfft.m 6021 2014-05-27 16:41:37Z vladimir $


%-This part if for creating a config branch that plugs into spm_cfg_eeg_tf
% Any parameters can be specified and they are then passed to the plugin
% when it's called.
%--------------------------------------------------------------------------
if nargin == 0
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
    
    mtmfft = cfg_branch;
    mtmfft.tag = 'mtmfft';
    mtmfft.name = 'Spectrum';
    mtmfft.val = {taper, freqres};
    
    res = mtmfft;
    
    return
elseif nargin < 3
    error('Three input arguments are required');
end

%-Defaults
%--------------------------------------------------------------------------
if ~isfield(S, 'taper')
    S.taper = 'sine';
end


dt = time(end) - time(1) + diff(time(1:2));


if ~isfield(S, 'frequencies') || isempty(S.frequencies)
    S.frequencies = (1/dt):max(1/dt, floor(dt)/dt):48;
end

if ~isfield(S, 'freqres')
    S.freqres = max(1/dt, floor(dt)/dt);
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

[spectrum,ntaper,freqoi] = ft_specest_mtmfft(data, time, 'taper', S.taper,'freqoi', S.frequencies,...
    'tapsmofrq', freqres, 'pad', pad, 'padtype', 'mirror', 'verbose', 0);

% This in principle should not happen
if length(unique(diff(freqoi)))>1 && max(abs(diff(unique(diff(freqoi)))))>1e-3
    error('Non-uniform output frequency axis')
end

res = [];
res.freq = linspace(freqoi(1), freqoi(end), length(freqoi));
res.time = [1 1]*mean(time);

if ndims(spectrum) == 3 && size(spectrum, 1)>1
    res.pow = spm_squeeze(nanmean(spectrum.*conj(spectrum), 1), 1);
elseif ndims(spectrum) == 3
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