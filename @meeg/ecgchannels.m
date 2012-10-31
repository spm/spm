function ind = ecgchannels(this)
% Return indices of ECG channels
% FORMAT ind = ecgchannels(this)
%
%  this      - MEEG object
%  ind       - row vector of indices of ECG channels
%
% See also eogchannels, emgchannels, meegchannels
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Christophe Phillips & Stefan Kiebel
% $Id: ecgchannels.m 5025 2012-10-31 14:44:13Z vladimir $

warning_flexible('ecgchannels method is deprecated. Use indchantype(D, ''ECG'')');
ind = indchantype(this, 'ECG');
