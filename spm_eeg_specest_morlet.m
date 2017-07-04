function res = spm_eeg_specest_morlet(S, data, time)
% Plugin for spm_eeg_tf implementing Morlet wavelet transform
% FORMAT res = spm_eeg_specest_morlet(S, data, time)
%
% S                     - input structure
% fields of S:
%    S.subsample   - factor by which to subsample the time axis (default - 1)
%  either
%    S.ncycles     - Morlet wavelet factor (default - 7)
%  or
%    S.timeres     - Fixed time window length in ms
%
%    S.frequencies - vector of frequencies (default - 0-48) at optimal frequency bins
%                            
% Output:
%  res - 
%   If no input is provided the plugin returns a cfg branch for itself
%
%   If input is provided:
%      res.fourier - the complex output of wavelet transform
%      res.time    - time axis
%      res.freq    - frequency axis
%__________________________________________________________________________
% Copyright (C) 2010-2017 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_specest_morlet.m 7129 2017-07-04 16:24:53Z guillaume $


%-This part if for creating a config branch that plugs into spm_cfg_eeg_tf
% Any parameters can be specified and they are then passed to the plugin
% when it's called.
%--------------------------------------------------------------------------
if nargin == 0
    subsample         = cfg_entry;
    subsample.tag     = 'subsample';
    subsample.name    = 'Subsample';
    subsample.strtype = 'n';
    subsample.num     = [1 1];
    subsample.val     = {1};
    subsample.help    = {'Set to N to subsample the time axis to every Nth sample (to reduce the dataset size).'};
    
    ncycles         = cfg_entry;
    ncycles.tag     = 'ncycles';
    ncycles.name    = 'Number of wavelet cycles';
    ncycles.strtype = 'n';
    ncycles.num     = [1 1];
    ncycles.val     = {7};
    ncycles.help    = {'Number of wavelet cycles (a.k.a. Morlet wavelet factor).',...
        'This parameter controls the time-frequency trade-off',...
        'Increasing it increases the frequency resolution at the expense of time resolution.'};
    
    timeres         = cfg_entry;
    timeres.tag     = 'timeres';
    timeres.name    = 'Fixed time window length';
    timeres.strtype = 'r';
    timeres.num     = [1 1];
    timeres.val     = {0};
    timeres.help    = {'Fixed time window for all frequencies.',...
        'Specify time window length in ms.',...
        'Default valued of 0 specifies variable time window length.'};
    
    morlet      = cfg_branch;
    morlet.tag  = 'morlet';
    morlet.name = 'Morlet wavelet transform';
    morlet.val  = {ncycles, timeres, subsample};
    
    res = morlet;
    
    return;
end

%-Defaults
%--------------------------------------------------------------------------
try, S.subsample; catch, S.subsample = 1; end
try, S.ncycles;   catch, S.ncycles = 7;   end

dt = time(end) - time(1);
if ~isfield(S, 'frequencies') || isempty(S.frequencies)
    S.frequencies = (1/dt):max(1/dt, floor(dt)/dt):48;
end

%-Generate wavelets
%--------------------------------------------------------------------------
if ~isfield(S, 'timeres') || S.timeres == 0
    M = spm_eeg_morlet(S.ncycles, 1000*diff(time(1:2)), S.frequencies);
else
    M = spm_eeg_morlet(S.ncycles, 1000*diff(time(1:2)), S.frequencies, 1000./S.timeres);
end

%-Data dimensions
%--------------------------------------------------------------------------
Nchannels    = size(data, 1);
Nsamples     = size(data, 2);
Nfrequencies = numel(M);

%-Initialize output struct
%--------------------------------------------------------------------------
res.freq    = S.frequencies;
res.time    = time(1:S.subsample:end);
res.fourier = zeros(Nchannels, Nfrequencies, length(res.time));

%-Compute wavelet transform
%--------------------------------------------------------------------------
for j = 1:Nchannels
    for i = 1:Nfrequencies
        tmp = conv(data(j, :), M{i});
        
        % time shift to remove delay
        tmp = tmp([1:Nsamples] + (length(M{i})-1)/2);
        
        tmp = tmp(1:S.subsample:end);
        
        res.fourier(j, i, :) = tmp;
    end
end

%% The following is often faster but requires more memory
%for i = 1:Nfrequencies
%    H = spm_convmtx(M{i}',Nsamples); % faster than conv for large Nchannels
%    tmp = data * H';
%    % subsample and time shift to remove delay
%    tmp = tmp(:,(1:S.subsample:Nsamples) + (length(M{i})-1)/2);
%    res.fourier(:,i,:) = tmp;
%end
