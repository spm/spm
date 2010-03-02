function res = spm_eeg_specest_hilbert(S, data, time)
% Plugin for spm_eeg_tf implementing spectral estimation using Hilbert transform
% FORMAT res = spm_eeg_specest_hilbert(S, data, time)
%
% S                     - input structure
% fields of S:
%    S.subsample   - factor by which to subsample the time axis (default - 1)
%    S.freqres     - frequency resolutions (plus-minus for each frequency, can
%                    be a vector with a value per frequency)
%    S.frequencies - vector of frequencies
%    S.order       - butterworth filter order (can be a vector with a value
%                    per frequency)
%                         
% Output:
%  res - 
%   If no input is provided the plugin returns a cfg branch for itself
%
%   If input is provided:
%      res.fourier - the complex output of wavelet transform
%      res.time    - time axis
%      res.freq    - frequency axis
%______________________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak based on the code contributed by Krish Singh
% $Id: spm_eeg_specest_hilbert.m 3742 2010-03-02 15:15:43Z vladimir $


%-This part if for creating a config branch that plugs into spm_cfg_eeg_tf
% Any parameters can be specified and they are then passed to the plugin
% when it's called.
%--------------------------------------------------------------------------
if nargin == 0 
    subsample = cfg_entry;
    subsample.tag = 'subsample';
    subsample.name = 'Subsample';
    subsample.strtype = 'n';
    subsample.num = [1 1];
    subsample.val = {1};
    subsample.help = {'Set to N to subsample the time axis to every Nth sample (to reduce the dataset size).'};
    
    freqres = cfg_entry;
    freqres.tag = 'freqres';
    freqres.name = 'Frequency resolution';
    freqres.strtype = 'r';
    freqres.num = [1 Inf];
    freqres.val = {0.5};
    freqres.help = {'Frequency resolution.',...
        'Note: 1 Hz resolution means plus-minus 1 Hz, i.e. 2 Hz badwidth',...
        'Either a single value or a vector of the same length as frequencies can be input'};
    
    order = cfg_entry;
    order.tag = 'order';
    order.name = 'Filter order';
    order.strtype = 'n';
    order.num = [1 Inf];
    order.val = {3};
    order.help = {'Butterworth filter order',...
        'Either a single value or a vector of the same length as frequencies can be input'};
    
    hilbert = cfg_branch;
    hilbert.tag = 'hilbert';
    hilbert.name = 'Hilbert transform';
    hilbert.val = {freqres, order, subsample};
    
    res = hilbert;
    
    return
elseif nargin < 3
    error('Three input arguments are required');
end

%-Defaults
%--------------------------------------------------------------------------
if ~isfield(S, 'subsample')
    S.subsample = 1;
end

dt = time(end) - time(1);

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
   
if ~isfield(S, 'order')
    S.order = 3;
end

if length(S.order) == 1
    order = S.order*ones(1, length(S.frequencies));
elseif length(S.order) == length(S.frequencies)
    order = S.order;
else
    error('Filter order should be either a scalar or a vector the same length as the number of frequencies.')
end

%-Data dimensions
%--------------------------------------------------------------------------
Nchannels = size(data, 1);
Nsamples = size(data, 2);
Nfrequencies = length(S.frequencies);

fsample = 1./diff(time(1:2));

%-Initialize output struct
%--------------------------------------------------------------------------
res = [];
res.freq = S.frequencies;
res.time = time(1:S.subsample:end);
res.fourier = zeros(Nchannels, Nfrequencies, length(res.time));

%-Compute wavelet transform
%--------------------------------------------------------------------------
for j = 1:Nchannels
    for i = 1:Nfrequencies
        highpass = S.frequencies(i) + freqres(i);
        lowpass  = S.frequencies(i) - freqres(i);
        
        if highpass < 1/dt
            error(sprintf('Frequency %f.2 Hz is too low to be estimated from data segment of %f.3 sec.', S.frequencies(i), dt));
        elseif lowpass < 1/dt
            padding = floor(fsample./highpass);
        else
            padding = floor(fsample./lowpass);
        end
            
        ind = [padding:-1:2 1:Nsamples (Nsamples-1):-1:(Nsamples-padding)];
        
        tmp = data(j, ind);
        
        if lowpass < 1/dt
            tmp = ft_preproc_lowpassfilter(tmp, fsample, highpass, order(i));
        else
            tmp = ft_preproc_bandpassfilter(tmp, fsample, [lowpass highpass], order(i));
        end
        
        tmp = spm_hilbert(tmp);
        
        % remove padding
        tmp = tmp(padding:(padding+Nsamples-1));
        
        tmp = tmp(1:S.subsample:end);
        
        res.fourier(j, i, :) = tmp;
    end
end
