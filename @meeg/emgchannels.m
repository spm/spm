function ind = emgchannels(this)
% Return indices of EMG channels
% FORMAT ind = emgchannels(this)
%
%  this      - MEEG object
%  ind       - row vector of indices of EMG channels
%
% See also eogchannels, ecgchannels, meegchannels
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Christophe Phillips & Stefan Kiebel
% $Id: emgchannels.m 5025 2012-10-31 14:44:13Z vladimir $

warning_flexible('emgchannels method is deprecated. Use indchantype(D, ''EMG'')');
ind = indchantype(this, 'EMG');
