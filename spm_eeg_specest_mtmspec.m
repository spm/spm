function res = spm_eeg_specest_mtmspec(S, data, time)
% Plugin for spm_eeg_tf using SPM implementation of multitaper method
% FORMAT res = spm_eeg_specest_mtmspec(S, data, time)
%
% S                     - input structure
% fields of S:
%    S.bandwidth   - time bandwidth parameter determining the degree of
%                    spectral smoothing (typically 3 or 4).
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
% $Id: spm_eeg_specest_mtmspec.m 4021 2010-07-28 12:43:16Z vladimir $


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
    
    bandwidth = cfg_entry;
    bandwidth.tag = 'bandwidth';
    bandwidth.name = 'Time bandwidth';
    bandwidth.strtype = 'n';
    bandwidth.num = [1 1];
    bandwidth.val = {3};
    bandwidth.help = {'Time bandwidth parameter (e.g. 3 or 4)'};
    
    mtmspec = cfg_branch;
    mtmspec.tag = 'mtmspec';
    mtmspec.name = 'SPM multitaper';
    mtmspec.val = {timeres, timestep, bandwidth};
    
    res = mtmspec;
    
    return
elseif nargin < 3
    error('Three input arguments are required');
end

%-Defaults
%--------------------------------------------------------------------------
if ~isfield(S, 'timeres')
    S.timeres = 400;
end

if ~isfield(S, 'timestep')
    S.timestep = 50;
end

if ~isfield(S, 'bandwidth')
    S.bandwidth = 3;
end

%-Data dimensions
%--------------------------------------------------------------------------
fsample = 1./diff(time(1:2));

timeres  = 1e-3*S.timeres;
timestep = 1e-3*S.timestep;

if timestep>timeres
    error('Time resolution should exceed time step');
end

%-Do the spectral analysis
%--------------------------------------------------------------------------
res = [];
[p, f, t] = spm_mmtspec(data', fsample, S.frequencies, timeres, timestep, S.bandwidth);

if size(data, 1) == 1
    res.pow = shiftdim(p, -1);
else
    res.pow = permute(p, [3 1 2]);
end

res.freq = f;

dt = (time(end)-time(1))./length(t);

res.time = (time(1)+dt/2):dt:time(end);

