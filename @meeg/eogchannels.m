function ind = eogchannels(this)
% Return indices of EOG channels
% FORMAT ind = eogchannels(this)
%
%  this      - MEEG object
%  ind       - row vector of indices of EOG channels
%
% See also ecgchannels, emgchannels, meegchannels
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: eogchannels.m 5025 2012-10-31 14:44:13Z vladimir $

warning_flexible('eogchannels method is deprecated. Use indchantype(D, ''EOG'')');
ind = indchantype(this, 'EOG');