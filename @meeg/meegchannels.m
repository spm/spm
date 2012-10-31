function ind = meegchannels(this, modality)
% Return indices of M/EEG channels
% FORMAT ind = meegchannels(this, modality)
%
%  this      - MEEG object
%  modality  - one of EEG, LFP, MEG (excluding planar), MEGPLANAR, MEEG [optional]
%
%  ind       - row vector of M/EEG channels
%
% See also eogchannels, ecgchannels, emgchannels
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: meegchannels.m 5025 2012-10-31 14:44:13Z vladimir $

warning_flexible('meegchannels method is deprecated. Use indchantype(D, modality)');
if nargin >1
    ind = indchantype(this, modality);
else
    ind = indchantype(this, 'MEEG');
end