function ind = emgchannels(this)
% Method for getting index vector of ecg channels
% FORMAT ind = emgchannels(this)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Christophe Phillips & Stefan Kiebel
% $Id: emgchannels.m 2668 2009-01-29 12:11:54Z christophe $

type = chantype(this);
ind = find(strcmpi('EMG', type));
ind = ind(:)'; % must be row to allow to use it as loop indices

