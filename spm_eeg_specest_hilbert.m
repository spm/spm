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
% $Id: spm_eeg_specest_hilbert.m 4463 2011-09-06 10:53:01Z vladimir $


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
    freqres.val = {0.5};
    freqres.help = {'Frequency resolution.',...
        'Note: 1 Hz resolution means plus-minus 1 Hz, i.e. 2 Hz badwidth',...
        'Either a single value or a vector of the same length as frequencies can be input'};        
     
    typ = cfg_menu;
    typ.tag = 'type';
    typ.name = 'Filter type';
    typ.labels = {'Butterworth', 'FIR'};
    typ.values = {'but', 'fir'};
    typ.val = {'but'};
    typ.help = {'Select the filter type.'};    
    
    dir = cfg_menu;
    dir.tag = 'dir';
    dir.name = 'Filter direction';
    dir.labels = {'Zero phase', 'Forward', 'Backward'};
    dir.values = {'twopass', 'onepass', 'onepass-reverse'};
    dir.val = {'twopass'};
    dir.help = {'Select the filter direction.'};
    
    order = cfg_entry;
    order.tag = 'order';
    order.name = 'Filter order';
    order.strtype = 'n';
    order.num = [1 Inf];
    order.val = {3};
    order.help = {'Butterworth filter order',...
        'Either a single value or a vector of the same length as frequencies can be input'};    
    
    flt = cfg_branch;
    flt.tag = 'filter';
    flt.name = 'Filter';
    flt.val = {typ dir order};
    
    polyorder = cfg_menu;
    polyorder.tag = 'polyorder';
    polyorder.name = 'Polynomial trend removal order';
    polyorder.labels = {'0', '1', '2', '3'};
    polyorder.values = {0, 1, 2, 3};
    polyorder.val = {1};
    polyorder.help = {'Order of polynome for trend removal (0 for no detrending)'};
    
    subsample = cfg_entry;
    subsample.tag = 'subsample';
    subsample.name = 'Subsample';
    subsample.strtype = 'n';
    subsample.num = [1 1];
    subsample.val = {1};
    subsample.help = {'Set to N to subsample the time axis to every Nth sample (to reduce the dataset size).'};
    
    hilbert = cfg_branch;
    hilbert.tag = 'hilbert';
    hilbert.name = 'Hilbert transform';
    hilbert.val = {freqres, flt, polyorder, subsample};
    
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

%-Data dimensions
%--------------------------------------------------------------------------
[res.fourier, res.freq, res.time] = ft_specest_hilbert(data, time, 'freqoi', S.frequencies, 'timeoi', time(1:S.subsample:end), 'width', S.freqres,...
    'filttype', S.filter.type, 'filtorder', S.filter.order, 'filtdir', S.filter.dir, 'polyremoval', S.polyorder, 'verbose', 0);